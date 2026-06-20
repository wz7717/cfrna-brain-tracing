# Bioinformatics Application Note change log

## Completed manuscript revisions

- Rewrote the manuscript as a final-candidate Bioinformatics Application Note draft centred on reusable Python and Streamlit software.
- Kept the core route unchanged: projected-VSD Network Top3 beam generation followed by logCPM-compatible resolution-group and exploratory exact-region reranking.
- Replaced the potentially post-hoc wording around route selection with development-analysis wording and fixed-route evaluation wording.
- Standardized terminology around `locked production route` and `default scoring route`; removed unnecessary use of competing labels such as formal or hybrid route except where describing component behaviour.
- Expanded the projected-VSD description to state that it is a projected variance-stabilized expression representation, uses fixed gene-wise slope/intercept/clipping parameters in `data/models/bo2023_reference_projector_linear_full.npz`, is single-sample compatible, is restricted to Network Top3 beam generation, and does not use target-cohort labels or cohort-level distributions.
- Tightened validation interpretation for AHBA, TCGA/BraTS and biofluid datasets to avoid claims of direct anatomical equivalence or clinical localization.
- Rewrote limitations to treat warnings as an intended resolution-aware output, not as software failure.
- Updated Supplementary information, Data availability, Figure 1 caption and alt text to match verified local files and the actual plotted panels.

## Verified repository contents

- Repository remote: `https://github.com/wz7717/cfrna-brain-tracing.git`
- README and installation instructions: `README.md`
- CLI quick-start command: `cfrna-tracing --help` from `setup.py`
- Streamlit quick-start command: `streamlit run streamlit_app.py`
- Command-line interface: `cli.py`
- Streamlit entry points: `streamlit_app.py`, `cfrna_tracing_app.py`, `app/`
- Network projection/scoring implementation: `core/network_tracing.py`
- Three-tier regional implementation: `core/bo2023_region_tracing.py`
- Projection utilities: `core/reference_projection.py`
- Unit tests: `tests/`
- Validation and benchmark scripts: `scripts/`, `benchmark_runner.py`
- Lightweight model artifacts: `data/models/`
- Supplementary material: `Bioinformatics_Application_Note_Supplementary_File_submission_20260619.pdf`, `Bioinformatics_Application_Note_Supplement_submission_20260619.md`
- Supplementary tables: `tables/TableS1_internal_validation_design.csv` through `tables/TableS6_claim_boundaries.csv`
- Figure 1 artwork: `figures/Figure1_cfRNA_BrainTrace_Bioinformatics_lowres.png` and `figures/Figure1_cfRNA_BrainTrace_Bioinformatics_highres_178mm.tif`

## Missing or unverified repository contents

- Added: root-level MIT `LICENSE` file for the project with copyright holder `王震`.
- Planned: release tag `v0.1.0`; not yet created locally because final submission assets are still being completed.
- Missing/unverified: Zenodo or Figshare archival DOI for the manuscript release.
- Stable Streamlit web-service URL supplied by author: `https://brain-cfrna-tracing.streamlit.app/`.
- Added: synthetic public example input/output package under `submission_ready_assets/example_io/`.
- Missing: separate public archival DOI for processed evaluation tables and figure source data.
- Missing/unverified: a dedicated `source_data/` directory for Figure 1. Supplementary tables exist locally and support the plotted values, but public archival source data still need to be finalized.

## Unresolved submission placeholders

- Corresponding author email: `wangzhen@cibr.ac.cn`
- `[ZENODO DOI REQUIRED BEFORE SUBMISSION]`
- Stable web-service URL: `https://brain-cfrna-tracing.streamlit.app/`.
- Software license: MIT License.
- Public example input/output: synthetic example package under `submission_ready_assets/example_io/`.
- Funding statement: `This work received no specific funding.`
- `[DOI REQUIRED BEFORE SUBMISSION]`

Unresolved: processed evaluation tables and figure source data still require public archival DOI before submission.

## Validation numbers

No validation numbers were changed. The final-candidate draft retains the supplied values:

- Projected-VSD Network LOSO Top1/Top3: 58.00%/91.58%
- Projected-VSD Network LOMO Top1/Top3: 53.72%/91.33%
- Locked production route LOSO Network Top3: 92.38%
- Locked production route LOMO Network Top3: 91.21%
- Resolution-group Top3 LOSO/LOMO: 72.36%/69.09%
- Exact-region Top3 LOSO/LOMO: 45.33%/42.36%
- AHBA Network Top1/Top3: 74.68%/94.42%
- AHBA resolution-group Top1/Top3: 36.26%/67.03%
- AHBA exact Top1/Top3: 24.18%/42.86%
- TCGA/BraTS Network Top3: 40.00%
- TCGA/BraTS broad-anatomy Top3: 64.62%

## Claim consistency findings

- The previous manuscript claim that the repository contains README, installation instructions, CLI, Streamlit entry points, tests, model artifacts and validation/benchmark scripts is supported by local files.
- The previous manuscript phrasing that implied an archived manuscript release was not supported by verified local repository contents. The final-candidate draft now keeps an explicit placeholder for DOI archival; the live Streamlit URL and synthetic example input/output package have been supplied.
- Figure 1 is not a placeholder. It contains a workflow panel and a validation bar-chart panel. It does not include a biofluid accuracy panel, so the caption and alt text now state that biofluid datasets lacking anatomical truth are stress tests and are not plotted as localization accuracy.
- README validation numbers are older/different from the final manuscript numbers, so the manuscript does not cite the README as the source for final validation values.

## Final Pre-Submission Checklist

| Item | Status | Note |
|---|---|---|
| Corresponding author email | Done | `wangzhen@cibr.ac.cn` supplied by author. |
| Funding statement | Done | `This work received no specific funding.` |
| License | Done | Root-level MIT `LICENSE` added; manuscript states MIT License. |
| Zenodo DOI | Missing | Archive the manuscript release and replace `[ZENODO DOI REQUIRED BEFORE SUBMISSION]`. |
| Stable web service URL | Done | `https://brain-cfrna-tracing.streamlit.app/` supplied by author. |
| Public example input/output | Done | Synthetic input/output package prepared under `submission_ready_assets/example_io/`. |
| Supplementary information file | Done | Supplementary PDF, Markdown and Tables S1-S6 exist locally. |
| Processed evaluation tables | Done locally | Tables S1-S6 exist locally; public DOI archival remains missing. |
| Figure source data | Needs author input | Local tables support plotted values, but public source-data archival remains missing. |
| Tagged release | Needs final action | Planned tag: `v0.1.0`; create after final assets are committed. |
| Repository README | Done | `README.md` exists. |
| Installation instructions | Done | Present in `README.md`. |
| CLI quick-start command | Done | `cfrna-tracing --help` is documented. |
| Streamlit quick-start instructions | Done | `streamlit run streamlit_app.py` is documented. |
| Unit tests | Done | `tests/` contains project tests. |
| Benchmark scripts | Done | `scripts/` and `benchmark_runner.py` exist. |
| Figure 1 finalized | Done locally | PNG and high-resolution TIF exist; final journal layout check still recommended. |
| Alt text finalized | Done | Updated to match actual Figure 1 content. |
| Data availability finalized | Missing | Public archival DOI and final repository release details remain unresolved. |
