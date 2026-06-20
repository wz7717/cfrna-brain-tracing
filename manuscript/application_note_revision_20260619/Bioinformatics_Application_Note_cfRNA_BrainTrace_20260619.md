# cfRNA-BrainTrace: hierarchical brain-origin inference from RNA-seq with a primate transcriptomic atlas

**Article type:** Application Note  
**Category:** Gene expression  
**Authors and affiliations:** [TO BE INSERTED BEFORE SUBMISSION]  
**Corresponding author:** [NAME AND EMAIL TO BE INSERTED BEFORE SUBMISSION]

## Abstract

### Summary

cfRNA-BrainTrace is a Python and Streamlit application for hierarchical brain-origin candidate ranking from RNA expression profiles. The current locked route maps each query from logCPM- or logTPM-compatible expression into Bo2023-like variance-stabilized space only for 10-class Network candidate generation, then reranks resolution-group and exact-region candidates in local logCPM-compatible expression space. This design avoids treating VSD projection as a new atlas and exposes anatomical resolution as an explicit output. In Bo2023 internal validation, the formal three-tier hybrid route achieved Network Top3 accuracy of 0.9238 in leave-one-sample-out validation and 0.9121 in leave-one-monkey-out validation. Resolution-group Top3 accuracy was 0.7236 and 0.6909, respectively, whereas exact-region Top3 remained lower at 0.4533 and 0.4236. In AHBA cross-species mapped-label validation, the hybrid route retained high Network Top3 accuracy (0.9442) and improved exact Top3 accuracy (0.4286) over both logCPM baseline and projected-VSD-only routes. TCGA/BraTS glioma tissue RNA-seq with MRI labels supported only coarse anatomical consistency. The software therefore provides reproducible hierarchical candidate ranking, confidence diagnostics and claim boundaries for atlas-based RNA source tracing.

### Availability and Implementation

Python 3.11+; command-line and Streamlit interfaces. Source code and documentation are hosted in a private development repository at **https://github.com/wz7717/cfrna-brain-tracing** and require public release, archival DOI and licence confirmation before submission. A stable live application URL is **[TO BE INSERTED BEFORE SUBMISSION]**.

### Contact

**[CORRESPONDING AUTHOR EMAIL REQUIRED BEFORE SUBMISSION]**

### Supplementary information

Supplementary data are available online.

## 1 Introduction

Cell-free and tissue RNA profiles can retain tissue-of-origin information, but anatomical inference within the brain is constrained by shared transcriptional programmes, atlas resolution and sample-domain shifts. A single exact-region label can therefore be misleading when the evidence supports only a broader candidate set. This is especially relevant for translational studies in which coarse anatomical consistency may be reproducible, while exact-region localization remains underdetermined.

cfRNA-BrainTrace implements a hierarchical route using the Bo2023 macaque brain transcriptomic atlas and a versioned anatomical dictionary. The software reports Network, resolution-group and exact-region candidates with coverage, entropy, margin and scope diagnostics. The principal endpoint is Network Top3 candidate generation; resolution groups are the main region-level reporting granularity; exact-region results are interpreted as local candidate rankings rather than deterministic calls.

## 2 Software description

### 2.1 Inputs and preprocessing

The application accepts tabular gene-expression input with gene symbols and abundance values. Metadata columns are preserved in exports when supplied. Gene identifiers are normalized, audited against the reference panel and intersected with model genes. External count-like or TPM-like matrices are converted to logCPM/logTPM-compatible space before downstream scoring. Matched-gene counts, non-zero coverage and out-of-scope warnings are reported for each sample.

### 2.2 Reference projection and hierarchical inference

The training reference remains the original Bo2023 expression resource: VSD/batch-removed expression for Network candidate generation and raw-feature-count-derived logCPM expression for local downstream reranking. A linear projector maps held-out or external query profiles into Bo2023-like VSD space for Network Top3 beam generation only. The projector fitted the Bo2023 training samples well at the sample level (global Pearson 0.9961; median sample Pearson 0.9969), but gene-level reconstruction was weaker, so the software does not claim precise VSD recovery for every gene.

After the projected-VSD Network Top3 beam is generated, candidate regions are restricted to those Networks. Resolution-group and exact-region candidates are then reranked using local discriminative genes in logCPM-compatible expression space. This hybrid route was selected because direct projected-VSD exact-region scoring was not consistently superior, especially in leave-one-monkey-out evaluation.

### 2.3 Interfaces and outputs

The Streamlit interface supports upload, metadata review, ranked Network and anatomical candidate display, downloadable tables and model warnings. The command-line interface supports scripted scoring, benchmark export and reproducible result bundles. Each report includes Top1 and Top3 candidates, rank scores, gene coverage, entropy, margins, route identifiers and resolution-specific labels. Accuracy is computed only when independently defined anatomical truth is available.

## 3 Implementation

cfRNA-BrainTrace is implemented in Python using NumPy, pandas and SciPy for scoring, scikit-learn utilities for evaluation, Plotly and Matplotlib for visualization, and Streamlit for the web interface. The Network projector route is implemented in `core/network_tracing.py`; the three-tier region route is implemented in `core/bo2023_region_tracing.py`; the web trigger and route text are exposed through `app/pages/tracing_page.py`. A smoke test on database sample `19R348` confirmed the production route `projected_vsd_network_top3_logcpm_resolution_local_exact`, Network overlap of 199/200 model genes and logCPM-based regional reranking.

## 4 Validation

### 4.1 Internal validation supports projected-VSD Network beam generation

In fold-local LOSO Network validation across 819 Bo2023 samples, projected VSD achieved Network Top1 accuracy of 0.5800 and Top3 accuracy of 0.9158, exceeding logCPM baseline (0.5556/0.8828) and native VSD (0.5311/0.8816). In strict LOMO validation, projected VSD achieved Network Top1 accuracy of 0.5372 and Top3 accuracy of 0.9133, again exceeding logCPM baseline and native VSD. These results support projected-VSD query representation for Network candidate generation.

### 4.2 Direct exact-region scoring is a diagnostic baseline, not the main endpoint

Direct exact-region projected-VSD scoring improved LOSO exact Top3 over logCPM baseline (0.3563 versus 0.2604), but in LOMO it remained below native VSD (0.3116 versus 0.3473). Exact-region Top1 and Top3 values were far below Network-level performance. These findings motivated the formal hybrid route and prevent claims that projected VSD is uniformly superior at exact-region resolution.

### 4.3 Formal three-tier hybrid validation

The complete three-tier route used projected-VSD Network Top3 beam generation followed by logCPM-compatible resolution-group and local exact-region reranking. In LOSO validation, Network, resolution-group and exact-region Top3 accuracies were 0.9238, 0.7236 and 0.4533, with median true ranks of 1.0, 2.0 and 4.0. In LOMO validation, the corresponding Top3 accuracies were 0.9121, 0.6909 and 0.4236. The gap between resolution-group and exact-region accuracy supports reporting resolution groups as the more defensible region-level endpoint.

### 4.4 External evaluations

In AHBA human brain RNA-seq with mapped labels, the hybrid route achieved Network Top1/Top3 accuracy of 0.7468/0.9442, group Top1/Top3 accuracy of 0.3626/0.6703 and exact Top1/Top3 accuracy of 0.2418/0.4286. Hybrid exact Top3 exceeded both logCPM baseline (0.3077) and projected-VSD-only scoring (0.2967), consistent with preserving Network-level projection gains while recovering local interpretability through logCPM reranking.

TCGA/BraTS glioma RNA-seq with corrected MRI-derived labels was used only for coarse anatomical consistency, because MRI truth labels are human atlas labels rather than Bo2023 macaque exact-region IDs. The hybrid route reached Network Top3 accuracy of 0.4000 and broad-anatomy Top3 accuracy of 0.6462, both above the logCPM baseline. GSE189919 projection established gene-overlap feasibility (15,622/21,668 projector genes; overlap fraction 0.7210) but was not used as a brain-origin accuracy validation because suitable anatomical truth labels were unavailable.

## 5 Conclusion

cfRNA-BrainTrace packages a three-tier atlas-based RNA source-tracing route into reproducible command-line and web interfaces. The evidence supports Network Top3 candidate generation and resolution-group reporting more strongly than exact-region Top1 localization. The software should therefore be used for hierarchical candidate ranking, quality control and transparent resolution auditing, not for direct clinical localization from unlabeled biofluid data.

## Funding

**[FUNDING INFORMATION TO BE INSERTED BEFORE SUBMISSION]**

## Conflict of Interest

None declared.

## Data availability

The Bo2023 macaque transcriptomic atlas, Allen Human Brain Atlas, TCGA/BraTS and GEO datasets are available from their original repositories under their original access conditions. Processed non-sensitive evaluation tables needed to reproduce manuscript figures should be included in the versioned software release.

## References

Bakas,S. *et al.* (2017) Advancing The Cancer Genome Atlas glioma MRI collections with expert segmentation labels and radiomic features. *Sci. Data*, **4**, 170117.

Bo,T. *et al.* (2023) Brain-wide and cell-specific transcriptomic insights into MRI-derived cortical morphology in macaque monkeys. *Nat. Commun.*, **14**, 1283.

Hawrylycz,M.J. *et al.* (2012) An anatomically comprehensive atlas of the adult human brain transcriptome. *Nature*, **489**, 391-399.

Vorperian,S.K. *et al.* (2022) Cell types of origin of the cell-free transcriptome. *Nat. Biotechnol.*, **40**, 855-861.

## Figure legend

**Fig. 1. Three-tier cfRNA-BrainTrace route and validation.** (A) Query expression is projected into Bo2023-like VSD space only for Network Top3 beam generation. (B) Resolution-group and exact-region candidates are reranked in local logCPM-compatible expression space. (C) Internal LOSO/LOMO validation supports Network Top3 and resolution-group reporting. (D) AHBA supports mapped-label external validation, whereas TCGA/BraTS supports only coarse anatomical consistency.

## Table

**Table 1. Principal route components and claim boundaries.**

| Component | Implementation | Manuscript claim boundary |
|---|---|---|
| Projector QC | Linear mapping from logCPM-compatible expression to Bo2023-like VSD space; sample-level global Pearson 0.9961 | Suitable for Network beam generation; not a new atlas and not gene-perfect VSD reconstruction |
| Network endpoint | Projected-VSD Top3 candidate beam | Main endpoint; LOSO Top3 0.9238 and LOMO Top3 0.9121 in formal route |
| Region-level endpoint | logCPM-compatible resolution-group reranking within Network beam | Primary region-level report; group Top3 exceeds exact Top3 in LOSO and LOMO |
| Exact-region output | Local logCPM exact-region reranking | Exploratory candidate list; do not present exact Top1 as deterministic localization |
| AHBA validation | Human brain RNA-seq with mapped labels | Cross-species mapped-label support; exact metrics only for evaluable labels |
| TCGA/BraTS validation | Glioma tissue RNA-seq plus MRI-derived human labels | Coarse anatomical consistency only; not Bo2023 exact-region validation |
