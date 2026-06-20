# cfRNA-BrainTrace: hierarchical brain-origin inference from RNA-seq with a primate transcriptomic atlas

**Article type:** Application Note  
**Category:** Gene expression  
**Authors and affiliations:** Author information to be supplied in the submission system  
**Corresponding author:** Corresponding author details to be supplied in the submission system

## Abstract

### Summary

cfRNA-BrainTrace is a Python and Streamlit application for hierarchical brain-origin candidate ranking from RNA expression profiles. Its formal route uses projected-VSD expression only to generate a broad Network Top3 beam, then returns to logCPM-compatible expression for resolution-group and exact-region reranking. The software reports anatomical resolution, confidence and scope diagnostics so that stable coarse signals are not overstated as deterministic exact-region calls.

### Availability and Implementation

Implemented in Python 3.11+ with command-line and Streamlit interfaces. Source code is hosted in the development repository **https://github.com/wz7717/cfrna-brain-tracing**. A versioned public release, archival DOI, licence statement and web service URL will accompany the final submission package.

### Contact

Corresponding author contact details will be provided in the submission system.

### Supplementary information

Supplementary data are available online.

## Text

### Motivation

RNA expression profiles can retain tissue-of-origin information, but within-brain source tracing is limited by transcriptional similarity among neighbouring regions, atlas granularity and domain shifts between reference tissue, tumour tissue and biofluid RNA. These constraints are especially important for atlas-based brain RNA tracing: a stable broad candidate may be recoverable even when exact-region localization is not. A single exact-region output can therefore imply more anatomical resolution than the data support. cfRNA-BrainTrace addresses this problem by returning a hierarchy of Network, resolution-group and exact-region candidates with explicit confidence, coverage and scope diagnostics.

### System and methods

The formal route was selected after comparing three behaviours observed during validation. First, projected-VSD query representation was consistently useful for broad Network candidate generation. Second, direct exact-region scoring in projected-VSD space was not consistently superior to native or logCPM-based alternatives. Third, local reranking within the Network beam improved interpretability at downstream anatomical resolutions. The resulting locked route is therefore a hybrid rather than a single-scale classifier.

Operationally, cfRNA-BrainTrace first normalizes the uploaded expression table into logCPM- or logTPM-compatible space and aligns genes to the reference panel. The query is then projected into Bo2023-like variance-stabilized space for 10-class Network scoring only. The top three Networks form a candidate beam; candidate regions outside this beam are excluded from downstream regional ranking. Within the retained beam, resolution groups are reranked using local logCPM-compatible expression, and exact regions are then reranked as lower-confidence local candidates. This route uses the Bo2023 macaque brain transcriptomic reference without converting the atlas into a new projected atlas. It is implemented in `core/network_tracing.py` and `core/bo2023_region_tracing.py`, with the same scoring core used by the command-line and Streamlit interfaces.

### Validation and outputs

Validation was organized around the same route logic used by the software. Network-level experiments first tested whether projected-VSD queries were suitable for generating the broad candidate beam. In fold-local leave-one-sample-out Network validation, projected VSD achieved Top1/Top3 accuracy of 58.00%/91.58%, exceeding logCPM baseline and native VSD routes. In strict leave-one-monkey-out validation, projected VSD reached 53.72%/91.33%. In contrast, direct exact-region projected-VSD scoring was not consistently superior, especially in leave-one-monkey-out evaluation, so exact-region scoring was not used as the sole endpoint.

The complete validation then tested the formal three-tier procedure end to end: projected-VSD Network Top3 beam generation, logCPM-compatible resolution-group reranking and logCPM-compatible exact-region reranking. This route achieved Network Top3 accuracy of 92.38% in leave-one-sample-out validation and 91.21% in leave-one-monkey-out validation. Resolution-group Top3 accuracy was 72.36% and 69.09%, respectively, whereas exact-region Top3 was 45.33% and 42.36%. The accuracy gradient across levels supports Network Top3 as the primary endpoint and resolution group as the most defensible region-level output. Exact-region results are retained as exploratory local candidate rankings rather than deterministic localization.

External analyses were matched to the resolution supported by their labels. In Allen Human Brain Atlas mapped-label validation, the hybrid route achieved Network Top1/Top3 accuracy of 74.68%/94.42%, resolution-group Top1/Top3 accuracy of 36.26%/67.03% and exact Top1/Top3 accuracy of 24.18%/42.86%; hybrid exact Top3 exceeded both logCPM baseline and projected-VSD-only scoring, supporting the use of projected VSD for the broad beam and logCPM expression for local reranking. In TCGA/BraTS glioma tissue RNA-seq with MRI-derived labels, results support only coarse anatomical consistency because MRI truth labels are human atlas labels rather than Bo2023 macaque exact-region identifiers. In that setting, hybrid Network Top3 was 40.00% and broad-anatomy Top3 was 64.62%. GSE189919 was used to verify projection feasibility rather than accuracy because anatomical truth labels were unavailable.

### Use and limitations

The software exports ranked candidates, matched-gene coverage, entropy, margin, route identifiers and resolution-specific warnings. Accuracy is calculated only when independent anatomical truth is available. Unlabelled biofluid predictions are reported as transfer stress tests or hypothesis-generating candidate rankings rather than clinical localization. The final release package is planned to include non-commercial access, test data, a stable public URL, versioned source code and an archival DOI.

## Funding

Funding information will be declared in the submission system.

## Conflict of Interest

None declared.

## Data availability

The Bo2023 macaque transcriptomic atlas, Allen Human Brain Atlas, TCGA/BraTS and GEO datasets are available from their original repositories under their original access conditions. Processed non-sensitive evaluation tables and figure source data are prepared for inclusion in the tagged public release and archive.

## References

Bakas,S. *et al.* (2017) Advancing The Cancer Genome Atlas glioma MRI collections with expert segmentation labels and radiomic features. *Sci. Data*, **4**, 170117.

Bo,T. *et al.* (2023) Brain-wide and cell-specific transcriptomic insights into MRI-derived cortical morphology in macaque monkeys. *Nat. Commun.*, **14**, 1283.

Hawrylycz,M.J. *et al.* (2012) An anatomically comprehensive atlas of the adult human brain transcriptome. *Nature*, **489**, 391-399.

Vorperian,S.K. *et al.* (2022) Cell types of origin of the cell-free transcriptome. *Nat. Biotechnol.*, **40**, 855-861.

## Figure legend

**Figure 1. cfRNA-BrainTrace three-tier route and validation evidence.** The query is projected into Bo2023-like VSD space only for Network Top3 beam generation; downstream resolution-group and exact-region reranking is performed in logCPM-compatible local expression space. Internal validation supports high Network Top3 accuracy in both LOSO and LOMO settings, while resolution-group accuracy is consistently higher than exact-region accuracy. AHBA supports mapped-label external validation, whereas TCGA/BraTS supports only coarse anatomical consistency.

**Alt text:** Multi-panel figure summarizing the cfRNA-BrainTrace workflow and validation. A workflow diagram shows query preprocessing, projected-VSD Network Top3 beam generation and logCPM local reranking. Bar charts show Network Top3 near 92% in internal validation, resolution-group Top3 around 69%-72% and lower exact-region Top3 around 42%-45%. External panels show AHBA Network Top3 of 94.42% and TCGA/BraTS broad-anatomy Top3 of 64.62%.
