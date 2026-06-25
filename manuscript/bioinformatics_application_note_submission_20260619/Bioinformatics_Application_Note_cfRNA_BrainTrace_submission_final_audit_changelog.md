# Final submission audit changelog

Date: 2026-06-25

## Manuscript changes

- Generated final clean Markdown and Word files from the current submission-ready manuscript text.
- The scientific mainline was retained. The formal LOSO Network denominator was corrected from the conditional 814-sample region-evaluable subset to all 819 Network-evaluable samples.
- Retained the three-tier route wording: uploaded RNA expression profiles are converted into logCPM/logTPM-compatible inputs, projected into Bo2023-like projected-VSD space only for 10-class macaque Network Top3 beam generation, and then reranked using logCPM-compatible local expression for resolution-group and exploratory exact-region candidates.
- Retained the AI-assisted editing disclosure.
- Retained the author placeholder `[AUTHOR NAMES AND AFFILIATIONS REQUIRED BEFORE SUBMISSION]` because complete author names and affiliations cannot be verified from the available files.

## README consistency

- README is synchronized with the manuscript route and validation numbers.
- README Status was updated to describe the v0.1.6 public submission release, research-use scope, non-clinical status, Zenodo archival and manuscript DOI citation.
- README no longer states that public repository URL, archived release DOI or OSI-approved license still need finalization.

## Current submission route

- Current route is consistently described as projected-VSD Network Top3 beam generation followed by logCPM-compatible resolution-group and exploratory exact-region reranking.
- Exact-region output is consistently described as exploratory, not deterministic localization.
- Biofluid datasets without anatomical truth are consistently treated as projection-feasibility or transfer stress tests, not localization-accuracy validation.

## Validation numbers

- Current submission numbers are consistent across the main manuscript, README, supplementary markdown, Tables S1-S6 and Figure 1 source data.
- Figure source data and Figure 1 labels use rounded display values consistent with the manuscript:
  - Network projected-VSD LOSO Top1/Top3: 58.00% / 91.58%.
  - Network projected-VSD LOMO Top1/Top3: 53.72% / 91.33%.
  - Locked three-tier route LOSO Network Top1/Top3: 58.24% / 92.19% (n=819).
  - Locked three-tier route LOMO Network Top3: 91.21%.
  - Resolution-group Top3 LOSO/LOMO: 72.36% / 69.09%.
  - Exact-region Top3 LOSO/LOMO: 45.33% / 42.36%.
  - AHBA mapped-label Network Top1/Top3: 74.68% / 94.42%.
  - AHBA resolution-group Top1/Top3: 36.26% / 67.03%.
  - AHBA exact Top1/Top3: 24.18% / 42.86%.
  - TCGA/BraTS Network Top3: 40.00%.
  - TCGA/BraTS broad-anatomy Top3: 64.62%.

## Legacy baseline

- README no longer reports legacy baseline validation numbers; it retains only a brief historical note that earlier baseline routes were used during development and are not part of the current submission route.
- The former formal LOSO Network Top3 value of 92.38% is explicitly labeled as a legacy denominator inconsistency.
- The alternate MRI-truth result of Network Top3 36.92% and broad Top3 80.00% is explicitly excluded from the locked Table S4 submission route.

## Zenodo DOI and archive

- Zenodo concept DOI: `https://doi.org/10.5281/zenodo.20773674`.
- Manuscript-associated v0.1.6 version DOI: `https://doi.org/10.5281/zenodo.20780280`.
- Manuscript and README identify both the fixed v0.1.6 version DOI and the project concept DOI.
- The Zenodo record is public and not draft/private.
- The Zenodo archive contains the GitHub source archive for the manuscript-associated release, including source code, README/release notes, synthetic example input/output, supplementary files, processed non-sensitive evaluation tables, figure source data, validation/benchmark scripts and license information.

## Streamlit demo

- Streamlit URL: `https://brain-cfrna-tracing.streamlit.app/`.
- Public HTTP reachability has been verified without login.
- Repository source inspection confirms upload/input-format guidance, example-compatible input handling, ranked candidate generation, diagnostics display and downloadable JSON/CSV outputs.
- The public demo text now explicitly explains marker coverage, entropy, score margin, scope warnings and the boundary that biofluid outputs without anatomical truth are transfer stress tests, not localization-validation results.

## Author block

- The corresponding author email `wangzhen@cibr.ac.cn` is retained.
- Full author names, departments, institution, city, postal code and country still require author input before submission.

## Supplementary files

- `Bioinformatics_Application_Note_Supplement_submission_20260619.md` exists and was expanded to include label harmonization, model artifact descriptions, projected-VSD projector description, example input/output description and locked evaluation timeline.
- `Bioinformatics_Application_Note_Supplementary_File_submission_20260619.pdf` was rebuilt from the updated supplementary markdown.
- Tables S1-S6 exist.
- Table S5 was corrected to index Table S5 and Table S6 explicitly.

## Tables S1-S6

- Table S1: internal validation design.
- Table S2: internal validation results.
- Table S3: external validation design and allowed conclusions.
- Table S4: external validation results.
- Table S5: figure/table index.
- Table S6: claim boundaries.
- Tables are consistent with the current manuscript route and validation numbers.

## Figure 1, caption and alt text

- Figure 1 contains query preprocessing, projected-VSD Network Top3 beam generation, logCPM-compatible reranking, internal LOSO/LOMO resolution gradient, AHBA mapped-label validation and TCGA/BraTS coarse-consistency results.
- Biofluid analyses are not plotted as localization accuracy.
- Caption and alt text are consistent with the figure and manuscript route.

## Figure source data

- Figure source data exist in `submission_ready_assets/source_data/`.
- `Figure1_validation_summary.csv` is the authoritative numeric source for Figure 1.
- `README_source_data.md` identifies both the manuscript-associated version DOI and the project concept DOI.
- Source-data tables are consistent with the manuscript validation numbers.

## Repository paths

All manuscript-mentioned paths were checked and exist:

- `core/network_tracing.py`
- `core/bo2023_region_tracing.py`
- `cli.py`
- `streamlit_app.py`
- `cfrna_tracing_app.py`
- `app/`
- `data/models/`
- `data/models/bo2023_reference_projector_linear_full.npz`
- `tests/`
- `scripts/`
- `benchmark_runner.py`
- `submission_ready_assets/example_io/`

## Projected-VSD implementation

- `data/models/bo2023_reference_projector_linear_full.npz` exists.
- The projector file contains `slope`, `intercept`, `clip_low` and `clip_high` arrays.
- Code inspection confirms single-sample projection through fixed reference-derived parameters.
- No target-cohort labels or target-cohort distribution fitting are used by the default projection path.
- Projected VSD is used for Network scoring/beam generation, while downstream resolution-group and exact-region reranking use logCPM-compatible local expression.

## Reproduction status

- Local unit tests: `17 passed`.
- Victor controlled-data workflow: completed successfully on 2026-06-25.
- Formal LOSO: Network `n=819`, Top1/Top3 58.24%/92.19%; resolution-group and exact-region `n=814`.
- Independent projected-VSD Network LOMO: `n=819`, Top1/Top3 53.72%/91.33%.
- Formal LOMO: Network `n=819`, Top1/Top3 57.75%/91.21%; resolution-group and exact-region `n=812`.
- AHBA mapped-label, TCGA/BraTS coarse-consistency and GSE189919 projection-feasibility results matched the submission tables.
- GSE189919's optional legacy algorithm audit remains a non-blocking technical limitation and does not support a localization-accuracy claim.
- The complete command, input, output, target, actual, consistency and difference-reason inventory is written to `output/v016_jupyter_validation/jupyter_validation_summary.csv`.

## AI disclosure and cover letter

- The manuscript retains `AI-assisted editing disclosure`.
- The cover letter should also disclose AI-assisted language editing, manuscript-format checking and code-documentation review support.

## Unresolved inconsistencies

- The formal LOSO denominator inconsistency was corrected: Network uses all
  819 samples, while resolution-group and exact-region metrics use the 814
  reference-supported folds. The former 92.38% Network Top3 value is legacy.
- Formal LOSO and formal LOMO use the same three-tier architecture but different
  fold-local Network construction; this is now disclosed rather than described
  as algorithmically identical.
- The only unresolved pre-submission item is the complete author/affiliation block, which needs author input before submission.

## Final recommendation

Ready after author input.
