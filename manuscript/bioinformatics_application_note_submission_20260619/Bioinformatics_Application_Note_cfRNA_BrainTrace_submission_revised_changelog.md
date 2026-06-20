# Revised Application Note change log

**Superseded change-log notice:** This file is retained for provenance only. The current submission audit is `Bioinformatics_Application_Note_cfRNA_BrainTrace_submission_final_audit_changelog.md`.

## Main wording changes

- Retitled the manuscript to emphasize brain-origin candidate ranking from RNA expression profiles rather than validated clinical cfRNA localization.
- Rewrote the Summary to state the locked route, main validation results and clinical/localization boundary.
- Replaced future-tense availability claims with verified repository contents and, at that stage, temporary missing-item markers that have since been resolved in the current final submission files.
- Defined Network as the 10-class macaque functional-anatomical source space and projected-VSD as projected variance-stabilized expression representation.
- Clarified that the reference-fitted linear projector uses stored slope/intercept/clipping parameters and does not use target-cohort labels or target-cohort distribution information in the default route.
- Added a separate Implementation section covering `core/network_tracing.py`, `core/bo2023_region_tracing.py`, `cli.py`, Streamlit entry points, model artifacts, tests and validation/benchmark scripts.
- Reframed validation around label-supported resolution: internal Network beam, formal three-tier route, AHBA mapped-label validation, TCGA/BraTS coarse consistency and GSE189919/biofluid transfer stress tests.
- Strengthened Use and limitations, including sparse profiles, marker coverage, entropy, margin, out-of-scope anatomy and domain shift warnings.

## Previously unresolved items, now resolved in the current final submission files

- Corresponding author email: resolved as `wangzhen@cibr.ac.cn`.
- Zenodo DOI: resolved as `https://doi.org/10.5281/zenodo.20773674`.
- Stable web service URL: resolved as `https://brain-cfrna-tracing.streamlit.app/`.
- License: resolved as MIT License.
- Example input/output: resolved under `submission_ready_assets/example_io/`.
- Funding: resolved as no specific funding.
- Processed table and figure-source archival DOI: resolved through the Zenodo concept DOI above.

## Repository items found

- README documentation and installation/run instructions: `README.md`
- Streamlit entry points: `streamlit_app.py`, `cfrna_tracing_app.py`, `app/`
- Command-line interface: `cli.py`; package entry point `cfrna-tracing=cli:main` in `setup.py`
- Tests: `tests/`
- Lightweight model artifacts: `data/models/`
- Validation and benchmark scripts: `scripts/`, `benchmark_runner.py`

## Repository items not found or not verified

- Historical note only; these items are resolved in the current final submission audit.

## Validation numbers

No validation numbers were changed. The revised draft retains the provided values for internal projected-VSD Network validation, formal three-tier LOSO/LOMO validation, AHBA mapped-label validation and TCGA/BraTS coarse anatomical consistency.
