# Revised Application Note change log

## Main wording changes

- Retitled the manuscript to emphasize brain-origin candidate ranking from RNA expression profiles rather than validated clinical cfRNA localization.
- Rewrote the Summary to state the locked route, main validation results and clinical/localization boundary.
- Replaced future-tense availability claims with verified repository contents plus explicit placeholders for DOI, web service URL, license, example input/output and corresponding author email.
- Defined Network as the 10-class macaque functional-anatomical source space and projected-VSD as projected variance-stabilized expression representation.
- Clarified that the reference-fitted linear projector uses stored slope/intercept/clipping parameters and does not use target-cohort labels or target-cohort distribution information in the default route.
- Added a separate Implementation section covering `core/network_tracing.py`, `core/bo2023_region_tracing.py`, `cli.py`, Streamlit entry points, model artifacts, tests and validation/benchmark scripts.
- Reframed validation around label-supported resolution: internal Network beam, formal three-tier route, AHBA mapped-label validation, TCGA/BraTS coarse consistency and GSE189919/biofluid transfer stress tests.
- Strengthened Use and limitations, including sparse profiles, marker coverage, entropy, margin, out-of-scope anatomy and domain shift warnings.

## Remaining placeholders

- `[CORRESPONDING AUTHOR EMAIL REQUIRED BEFORE SUBMISSION]`
- `[ZENODO DOI REQUIRED BEFORE SUBMISSION]`
- `[STABLE WEB SERVICE URL REQUIRED BEFORE SUBMISSION]`
- `[LICENSE REQUIRED BEFORE SUBMISSION]`
- `[EXAMPLE INPUT/OUTPUT REQUIRED BEFORE SUBMISSION]`
- `[FUNDING INFORMATION REQUIRED BEFORE SUBMISSION]`
- `[DOI REQUIRED BEFORE SUBMISSION]`

## Repository items found

- README documentation and installation/run instructions: `README.md`
- Streamlit entry points: `streamlit_app.py`, `cfrna_tracing_app.py`, `app/`
- Command-line interface: `cli.py`; package entry point `cfrna-tracing=cli:main` in `setup.py`
- Tests: `tests/`
- Lightweight model artifacts: `data/models/`
- Validation and benchmark scripts: `scripts/`, `benchmark_runner.py`

## Repository items not found or not verified

- Root-level OSI-approved license file
- Git tag or tagged release
- Zenodo/Figshare archival DOI
- Confirmed live Streamlit/web service URL
- Public example input/output files suitable for submission

## Validation numbers

No validation numbers were changed. The revised draft retains the provided values for internal projected-VSD Network validation, formal three-tier LOSO/LOMO validation, AHBA mapped-label validation and TCGA/BraTS coarse anatomical consistency.
