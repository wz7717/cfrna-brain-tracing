# Final pre-submission check changelog

Date: 2026-06-21

## Manuscript revisions

- Replaced the generic title-page author placeholder with `[AUTHOR NAMES AND AFFILIATIONS REQUIRED BEFORE SUBMISSION]` while retaining the confirmed corresponding email `wangzhen@cibr.ac.cn`.
- Compressed the Summary to a more Application Note-like software summary while preserving the locked route, main validation numbers and claim boundary.
- Added an `AI-assisted editing disclosure` section after `Conflict of Interest`.
- Preserved the restrained route-development wording to avoid implying post-hoc independent confirmation.
- Confirmed that the manuscript describes projected VSD as restricted to Network Top3 beam generation and logCPM-compatible expression as the basis for downstream resolution-group and exploratory exact-region reranking.

## README consistency

- Updated README workflow from the legacy `Bo2023 VSD reference -> fold-selected 200 genes -> Pearson correlation -> Top-3 pairwise rescue` description to the current submission route.
- Updated README validation summary to the current submission numbers.
- Moved the old 55.8%/88.0%, 53.2%/86.7%, 32.6%/55.4% and 75.4%/83.1% values into a `Legacy baseline / previous route` subsection.
- Updated README status to state that the repository is a research software release under the MIT License with Zenodo archival.

## Supplementary files

- Confirmed the supplementary markdown exists.
- Added `Model development and locked evaluation timeline` to the supplementary methods.
- Updated the supplementary figure/table description from Tables S1-S5 to Tables S1-S6.
- Rebuilt `Bioinformatics_Application_Note_Supplementary_File_submission_20260619.pdf` from the updated supplementary markdown.
- Confirmed local Tables S1-S6 exist and include validation design, internal validation results, external validation design/results, figure/table index and claim boundaries.
- PDF text extraction confirmed the rebuilt supplement includes the locked evaluation timeline. Page-image rendering was not performed because `pdftoppm` was not available in the current environment.

## Zenodo DOI

- DOI checked: `https://doi.org/10.5281/zenodo.20773674`.
- Version record checked: `https://doi.org/10.5281/zenodo.20775154`.
- Public accessibility was verified without login for the v0.1.2 record on 2026-06-21.
- The Zenodo API record reports title `wz7717/cfrna-brain-tracing: cfRNA-BrainTrace v0.1.2`, DOI `10.5281/zenodo.20775154`, concept DOI `10.5281/zenodo.20773674` and source archive `wz7717/cfrna-brain-tracing-v0.1.2.zip`.
- The manuscript retains the Zenodo concept DOI so the citation resolves to the latest archived release.

## Streamlit demo

- Demo URL checked: `https://brain-cfrna-tracing.streamlit.app/`.
- Public reachability was verified by HTTP 200 response on 2026-06-21.
- The manuscript retains the live demo URL.

## License

- `LICENSE` exists in the repository.
- README and manuscript state the MIT License.

## GitHub release/tag

- GitHub repository: `https://github.com/wz7717/cfrna-brain-tracing`.
- Submission release/tag: `v0.1.2`.
- Release assets include `cfRNA-BrainTrace_v0.1.2_source_data.zip` and `cfRNA-BrainTrace_v0.1.2_example_io.zip`.
- GitHub API access to the v0.1.2 release was verified on 2026-06-21.

## Example input/output

- Example files exist under `submission_ready_assets/example_io/`.
- Release asset `cfRNA-BrainTrace_v0.1.2_example_io.zip` contains synthetic input templates, full 200-gene input, ranked candidate CSVs, JSON outputs and the generation script.

## Tests and benchmark scripts

- `tests/` exists and contains unit tests for Network scoring, VSD adaptation, marker-route behaviour, region-resolution annotations and upload metadata.
- `benchmark_runner.py` exists.
- Validation/export scripts exist under `scripts/`.

## Figure source data

- Figure source data exist under `submission_ready_assets/source_data/` and in `cfRNA-BrainTrace_v0.1.2_source_data.zip`.
- Figure 1 low-resolution PNG was visually inspected and is consistent with the caption: projected VSD is restricted to Network beam generation, downstream reranking is logCPM-compatible, and TCGA/BraTS is shown as coarse consistency rather than exact-region accuracy.

## Validation-number consistency

- No current-route validation numbers were changed in the manuscript.
- README was corrected because it still presented legacy baseline values as if they were current-route results.
- Current submission numbers retained:
  - Network projected-VSD LOSO Top1/Top3: 58.00% / 91.58%.
  - Network projected-VSD LOMO Top1/Top3: 53.72% / 91.33%.
  - Locked three-tier LOSO Network Top3: 92.38%.
  - Locked three-tier LOMO Network Top3: 91.21%.
  - Resolution-group Top3 LOSO/LOMO: 72.36% / 69.09%.
  - Exact-region Top3 LOSO/LOMO: 45.33% / 42.36%.
  - AHBA Network Top1/Top3: 74.68% / 94.42%.
  - AHBA resolution-group Top1/Top3: 36.26% / 67.03%.
  - AHBA exact-region Top1/Top3: 24.18% / 42.86%.
  - TCGA/BraTS Network Top3: 40.00%.
  - TCGA/BraTS broad-anatomy Top3: 64.62%.

## Unresolved pre-submission items

- Full author names, departments, institution, city, postal code and country must still be supplied by the authors.
- Cover letter should disclose AI-assisted language editing, manuscript-format checking and code-documentation review support.
- The updated README and manuscript-ready checked files should be pushed and, if desired, archived as a new GitHub/Zenodo release after author approval.

## Final recommendation

Ready after author input.
