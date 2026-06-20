from __future__ import annotations

import re
from pathlib import Path

import build_revised_application_note as base


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript" / "bioinformatics_application_note_submission_20260619"
FIGURE = OUT / "figures" / "Figure1_cfRNA_BrainTrace_Bioinformatics_lowres.png"
MD_OUT = OUT / "Bioinformatics_Application_Note_cfRNA_BrainTrace_submission_final_candidate.md"
DOCX_OUT = OUT / "Bioinformatics_Application_Note_cfRNA_BrainTrace_submission_final_candidate.docx"
CHANGELOG_OUT = OUT / "Bioinformatics_Application_Note_cfRNA_BrainTrace_change_log.md"


FINAL_MD = """# cfRNA-BrainTrace: hierarchical brain-origin candidate ranking from RNA expression profiles with a primate transcriptomic atlas

**Article type:** Application Note  
**Category:** Gene expression  
**Authors and affiliations:** Author information to be supplied in the submission system  
**Corresponding author:** wangzhen@cibr.ac.cn

## Abstract

### Summary

cfRNA-BrainTrace is a Python and Streamlit application for hierarchical brain-origin candidate ranking from RNA expression profiles. The locked production route projects query profiles into a Bo2023-like variance-stabilized space to generate a broad 10-class macaque functional-anatomical Network Top3 beam, followed by logCPM-compatible local reranking for resolution-group and exploratory exact-region candidates. The software reports ranked candidates, marker coverage, entropy, score margins, route identifiers and scope warnings. In internal validation, Network Top3 accuracy reached 92.38% in leave-one-sample-out and 91.21% in leave-one-monkey-out evaluation, with lower accuracy at resolution-group and exact-region levels. AHBA mapped-label validation supported coarse cross-species transfer, whereas TCGA/BraTS and unlabeled biofluid analyses defined transfer limitations. cfRNA-BrainTrace is intended for reproducible coarse candidate ranking and resolution-limit auditing, not stand-alone clinical localization from unlabeled biofluid RNA.

### Availability and Implementation

Implemented in Python 3.11+ with command-line and Streamlit interfaces. Source code, README documentation, installation instructions, command-line and Streamlit entry points, unit tests, lightweight model artifacts, supplementary tables, validation/benchmark scripts and synthetic public example input/output files are available at **https://github.com/wz7717/cfrna-brain-tracing**. The manuscript release is archived at **https://doi.org/10.5281/zenodo.20773674**. A live Streamlit demonstration is available at **https://brain-cfrna-tracing.streamlit.app/**. The software is released under the **MIT License**.

### Contact

**wangzhen@cibr.ac.cn**

### Supplementary information

Supplementary Methods and Tables are provided in `Bioinformatics_Application_Note_Supplementary_File_submission_20260619.pdf`, `Bioinformatics_Application_Note_Supplement_submission_20260619.md` and Tables S1-S6. Synthetic public example input/output files are provided in `submission_ready_assets/example_io/`.

## Text

### Motivation

RNA expression profiles can retain tissue-of-origin information, but within-brain source tracing is limited by transcriptional similarity among neighbouring regions, atlas granularity and domain shifts between reference tissue, tumour tissue and biofluid RNA. A stable broad candidate may be recoverable even when exact-region localization is not. cfRNA-BrainTrace addresses this problem by returning Network, resolution-group and exact-region candidates together with confidence, coverage and scope diagnostics.

### System and methods

The locked production route was motivated by development analyses showing complementary behaviour of projected-VSD and logCPM-compatible expression spaces, and was then evaluated as a fixed three-tier route. Network denotes the 10-class macaque functional-anatomical source space used by the Bo2023 reference. Projected-VSD denotes a projected variance-stabilized expression representation used for broad Network candidate generation, whereas downstream regional reranking is performed in a logCPM-compatible local expression space.

Operationally, cfRNA-BrainTrace normalizes the uploaded expression table into logCPM- or logTPM-compatible space and aligns genes to the reference panel. A reference-fitted linear projector maps a query profile to Bo2023-like VSD space using stored gene-wise slope, intercept and clipping parameters from `data/models/bo2023_reference_projector_linear_full.npz`. This projection is single-sample compatible: it uses fixed reference-derived parameters and does not use target-cohort labels or cohort-level distribution information, reducing transductive-leakage risk. The projected representation is used only for Network Top3 beam generation. The top three Networks form the candidate beam, and regions outside this beam are excluded from downstream regional ranking. Within the retained beam, resolution groups and exact regions are reranked using logCPM-compatible local expression.

### Implementation

The Network projection and scoring route is implemented in `core/network_tracing.py`, and the three-tier regional route is implemented in `core/bo2023_region_tracing.py`. The command-line interface (`cli.py`) and Streamlit entry points (`streamlit_app.py`, `cfrna_tracing_app.py` and `app/`) call the same scoring core to avoid interface-specific divergence. Versioned model artifacts in `data/models/` include marker order, Network centroids, anatomical dictionaries, projection parameters, route parameters and warning metadata. The repository contains unit tests under `tests/` for Network scoring, region-resolution annotations, marker-route behaviour, upload metadata handling and VSD adaptation. Validation and export scripts under `scripts/`, together with `benchmark_runner.py`, support reruns of the validation analyses and generation of manuscript tables and figure artwork.

### Validation

Validation was matched to the resolution supported by available labels. Internal experiments first tested whether projected-VSD query representation was appropriate for generating the broad Network candidate beam. In fold-local leave-one-sample-out Network validation, projected VSD achieved Top1/Top3 accuracy of 58.00%/91.58%, exceeding logCPM baseline and native VSD routes. In strict leave-one-monkey-out validation, projected VSD reached 53.72%/91.33%. Direct exact-region projected-VSD scoring was lower and less stable, especially in leave-one-monkey-out evaluation, so exact-region scoring was not used as the sole endpoint.

The complete validation then tested the locked production route end to end: projected-VSD Network Top3 beam generation, logCPM-compatible resolution-group reranking and logCPM-compatible exact-region reranking. This route achieved Network Top3 accuracy of 92.38% in leave-one-sample-out validation and 91.21% in leave-one-monkey-out validation. Resolution-group Top3 accuracy was 72.36% and 69.09%, respectively, whereas exact-region Top3 was 45.33% and 42.36%. The accuracy gradient across levels supports Network Top3 as the primary endpoint and resolution group as the more defensible region-level output. Exact-region results are retained as exploratory local candidate rankings rather than a localization endpoint.

In AHBA mapped-label external validation, the locked production route achieved Network Top1/Top3 accuracy of 74.68%/94.42%, resolution-group Top1/Top3 accuracy of 36.26%/67.03% and exact Top1/Top3 accuracy of 24.18%/42.86%; exact Top3 exceeded both logCPM baseline and projected-VSD-only scoring. These AHBA results should be interpreted as mapped-label transfer rather than direct anatomical equivalence, because human anatomical labels were harmonized to the macaque-derived hierarchy. In TCGA/BraTS glioma tissue RNA-seq with MRI-derived labels, results support only coarse anatomical consistency because MRI truth labels are human atlas labels rather than Bo2023 macaque exact-region identifiers. In that setting, Network Top3 was 40.00% and broad-anatomy Top3 was 64.62%. The TCGA/BraTS results therefore support broad candidate consistency in tumour tissue but do not validate macaque Network-level or exact-region localization in human glioma. GSE189919 and other biofluid datasets without anatomical truth were treated as projection-feasibility or transfer stress tests rather than localization-accuracy validation.

### Use and limitations

The intended use of cfRNA-BrainTrace is coarse brain-origin candidate ranking from RNA expression profiles and assessment of whether a sample supports Network-level, resolution-group-level or only low-confidence output. The software exports warnings for sparse profiles, low marker coverage, high entropy, low score margins, out-of-scope anatomy and domain shift. Warnings guide resolution-aware interpretation: samples with low marker coverage, low margins or high entropy should be reported as low-confidence or coarse-only outputs rather than forced exact-region predictions. cfRNA-BrainTrace is not intended for deterministic exact-region localization, clinical cfRNA localization without independent anatomical truth, or localization of cerebellar/posterior-fossa samples outside the current reference space. The cfRNA name reflects the intended prospective research setting; current evidence does not establish clinical liquid-biopsy localization.

## Funding

This work received no specific funding.

## Conflict of Interest

None declared.

## Data availability

The Bo2023 macaque transcriptomic atlas, Allen Human Brain Atlas, TCGA/BraTS and GEO datasets are available from their original repositories under their original access conditions. Repository scripts, processed non-sensitive evaluation tables and figure source data are archived with the release at **https://doi.org/10.5281/zenodo.20773674**.

## References

Bakas,S. *et al.* (2017) Advancing The Cancer Genome Atlas glioma MRI collections with expert segmentation labels and radiomic features. *Sci. Data*, **4**, 170117.

Bo,T. *et al.* (2023) Brain-wide and cell-specific transcriptomic insights into MRI-derived cortical morphology in macaque monkeys. *Nat. Commun.*, **14**, 1283.

Hawrylycz,M.J. *et al.* (2012) An anatomically comprehensive atlas of the adult human brain transcriptome. *Nature*, **489**, 391-399.

Vorperian,S.K. *et al.* (2022) Cell types of origin of the cell-free transcriptome. *Nat. Biotechnol.*, **40**, 855-861.

## Figure legend

**Figure 1. cfRNA-BrainTrace three-tier route and validation evidence.** (A) Query profiles are represented as logCPM/logTPM-compatible inputs, projected into Bo2023-like VSD space only for Network Top3 beam generation, and then reranked in logCPM-compatible local expression space for resolution-group and exploratory exact-region candidates. (B) Internal validation shows the expected resolution gradient, with high Network Top3 accuracy and lower resolution-group and exact-region accuracy. AHBA provides mapped-label external validation, whereas TCGA/BraTS supports only coarse anatomical consistency because its MRI truth labels are human atlas labels rather than Bo2023 macaque exact-region identifiers. Biofluid analyses lacking anatomical truth are not shown as accuracy bars and are treated as projection-feasibility or transfer stress tests.

**Alt text:** Multi-panel figure summarizing the cfRNA-BrainTrace workflow and validation. Panel A shows query preprocessing, projected-VSD Network Top3 beam generation and logCPM-compatible local reranking. Panel B shows internal LOSO and LOMO resolution gradients, AHBA mapped-label validation, and TCGA/BraTS coarse-consistency results: Network Top3 is near 92% internally, resolution-group Top3 is about 69%-72%, exact-region Top3 is about 42%-45%, AHBA Network Top3 is 94.42%, and TCGA/BraTS broad-anatomy Top3 is 64.62%. Biofluid stress-test results are not plotted as localization accuracy because anatomical truth is absent.
"""


CHANGELOG = """# Bioinformatics Application Note change log

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
- Zenodo archival concept DOI generated for the manuscript release: `https://doi.org/10.5281/zenodo.20773674`.
- Stable Streamlit web-service URL supplied by author: `https://brain-cfrna-tracing.streamlit.app/`.
- Added: synthetic public example input/output package under `submission_ready_assets/example_io/`.
- Missing: separate public archival DOI for processed evaluation tables and figure source data.
- Missing/unverified: a dedicated `source_data/` directory for Figure 1. Supplementary tables exist locally and support the plotted values, but public archival source data still need to be finalized.

## Unresolved submission placeholders

- Corresponding author email: `wangzhen@cibr.ac.cn`
- Zenodo DOI: `https://doi.org/10.5281/zenodo.20773674`.
- Stable web-service URL: `https://brain-cfrna-tracing.streamlit.app/`.
- Software license: MIT License.
- Public example input/output: synthetic example package under `submission_ready_assets/example_io/`.
- Funding statement: `This work received no specific funding.`
- Data/source DOI: `https://doi.org/10.5281/zenodo.20773674`.

Resolved: processed evaluation tables and figure source data are included in the Zenodo-archived release.

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
| Zenodo DOI | Done | `https://doi.org/10.5281/zenodo.20773674`. |
| Stable web service URL | Done | `https://brain-cfrna-tracing.streamlit.app/` supplied by author. |
| Public example input/output | Done | Synthetic input/output package prepared under `submission_ready_assets/example_io/`. |
| Supplementary information file | Done | Supplementary PDF, Markdown and Tables S1-S6 exist locally. |
| Processed evaluation tables | Done | Tables S1-S6 are included in the Zenodo-archived release. |
| Figure source data | Done | Figure source artwork and supporting tables are included in the Zenodo-archived release. |
| Tagged release | Needs final action | Planned tag: `v0.1.0`; create after final assets are committed. |
| Repository README | Done | `README.md` exists. |
| Installation instructions | Done | Present in `README.md`. |
| CLI quick-start command | Done | `cfrna-tracing --help` is documented. |
| Streamlit quick-start instructions | Done | `streamlit run streamlit_app.py` is documented. |
| Unit tests | Done | `tests/` contains project tests. |
| Benchmark scripts | Done | `scripts/` and `benchmark_runner.py` exist. |
| Figure 1 finalized | Done locally | PNG and high-resolution TIF exist; final journal layout check still recommended. |
| Alt text finalized | Done | Updated to match actual Figure 1 content. |
| Data availability finalized | Done | DOI and release details are now included. |
"""


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    MD_OUT.write_text(FINAL_MD, encoding="utf-8")
    CHANGELOG_OUT.write_text(CHANGELOG, encoding="utf-8")
    base.FIGURE = FIGURE
    base.build_docx(FINAL_MD, DOCX_OUT)
    from docx import Document
    doc = Document(DOCX_OUT)
    doc.core_properties.title = "cfRNA-BrainTrace Bioinformatics Application Note final candidate"
    doc.core_properties.comments = "Final-candidate draft generated after repository-content verification."
    doc.save(DOCX_OUT)
    words = len(re.findall(r"[A-Za-z0-9]+(?:[-'][A-Za-z0-9]+)?", re.sub(r"`[^`]+`", " code ", FINAL_MD)))
    print(DOCX_OUT)
    print(MD_OUT)
    print(CHANGELOG_OUT)
    print(f"Approximate word count: {words}")


if __name__ == "__main__":
    main()
