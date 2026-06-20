# cfRNA-BrainTrace: hierarchical brain-origin candidate ranking from RNA expression profiles with a primate transcriptomic atlas

**Article type:** Application Note  
**Category:** Gene expression  
**Authors and affiliations:** Author information to be supplied in the submission system  
**Corresponding author:** wangzhen@cibr.ac.cn

## Abstract

### Summary

cfRNA-BrainTrace is a Python and Streamlit application for hierarchical brain-origin candidate ranking from RNA expression profiles. The locked production route projects query profiles into a Bo2023-like variance-stabilized space to generate a broad 10-class macaque functional-anatomical Network Top3 beam, followed by logCPM-compatible local reranking for resolution-group and exploratory exact-region candidates. The software reports ranked candidates, marker coverage, entropy, score margins, route identifiers and scope warnings. In internal validation, Network Top3 accuracy reached 92.38% in leave-one-sample-out and 91.21% in leave-one-monkey-out evaluation, with lower accuracy at resolution-group and exact-region levels. AHBA mapped-label validation supported coarse cross-species transfer, whereas TCGA/BraTS and unlabeled biofluid analyses defined transfer limitations. cfRNA-BrainTrace is intended for reproducible coarse candidate ranking and resolution-limit auditing, not stand-alone clinical localization from unlabeled biofluid RNA.

### Availability and Implementation

Implemented in Python 3.11+ with command-line and Streamlit interfaces. Source code, README documentation, installation instructions, command-line and Streamlit entry points, unit tests, lightweight model artifacts, supplementary tables, validation/benchmark scripts and synthetic public example input/output files are available at **https://github.com/wz7717/cfrna-brain-tracing**. The manuscript release is archived at **https://doi.org/10.5281/zenodo.20773675**. A live Streamlit demonstration is available at **https://brain-cfrna-tracing.streamlit.app/**. The software is released under the **MIT License**.

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

The Bo2023 macaque transcriptomic atlas, Allen Human Brain Atlas, TCGA/BraTS and GEO datasets are available from their original repositories under their original access conditions. Repository scripts, processed non-sensitive evaluation tables and figure source data are archived with the release at **https://doi.org/10.5281/zenodo.20773675**.

## References

Bakas,S. *et al.* (2017) Advancing The Cancer Genome Atlas glioma MRI collections with expert segmentation labels and radiomic features. *Sci. Data*, **4**, 170117.

Bo,T. *et al.* (2023) Brain-wide and cell-specific transcriptomic insights into MRI-derived cortical morphology in macaque monkeys. *Nat. Commun.*, **14**, 1283.

Hawrylycz,M.J. *et al.* (2012) An anatomically comprehensive atlas of the adult human brain transcriptome. *Nature*, **489**, 391-399.

Vorperian,S.K. *et al.* (2022) Cell types of origin of the cell-free transcriptome. *Nat. Biotechnol.*, **40**, 855-861.

## Figure legend

**Figure 1. cfRNA-BrainTrace three-tier route and validation evidence.** (A) Query profiles are represented as logCPM/logTPM-compatible inputs, projected into Bo2023-like VSD space only for Network Top3 beam generation, and then reranked in logCPM-compatible local expression space for resolution-group and exploratory exact-region candidates. (B) Internal validation shows the expected resolution gradient, with high Network Top3 accuracy and lower resolution-group and exact-region accuracy. AHBA provides mapped-label external validation, whereas TCGA/BraTS supports only coarse anatomical consistency because its MRI truth labels are human atlas labels rather than Bo2023 macaque exact-region identifiers. Biofluid analyses lacking anatomical truth are not shown as accuracy bars and are treated as projection-feasibility or transfer stress tests.

**Alt text:** Multi-panel figure summarizing the cfRNA-BrainTrace workflow and validation. Panel A shows query preprocessing, projected-VSD Network Top3 beam generation and logCPM-compatible local reranking. Panel B shows internal LOSO and LOMO resolution gradients, AHBA mapped-label validation, and TCGA/BraTS coarse-consistency results: Network Top3 is near 92% internally, resolution-group Top3 is about 69%-72%, exact-region Top3 is about 42%-45%, AHBA Network Top3 is 94.42%, and TCGA/BraTS broad-anatomy Top3 is 64.62%. Biofluid stress-test results are not plotted as localization accuracy because anatomical truth is absent.
