# cfRNA-BrainTrace: hierarchical brain-origin inference from RNA-seq with a primate transcriptomic atlas

**Article type:** Application Note  
**Category:** Gene expression  
**Authors and affiliations:** [TO BE INSERTED BEFORE SUBMISSION]  
**Corresponding author:** [NAME AND EMAIL TO BE INSERTED BEFORE SUBMISSION]

## Abstract

### Summary

cfRNA-BrainTrace is a Python and Streamlit application for hierarchical inference of brain-origin candidates from bulk or biofluid RNA expression profiles. The software aligns gene symbols to a primate brain reference, scores 10 broad anatomical Networks using a fold-selected 200-gene Pearson-correlation model, applies a low-margin pairwise rescue rule and reports Network, broad-anatomy, lobe and exploratory exact-region outputs with confidence diagnostics. Command-line and web interfaces generate ranked predictions, audit tables and publication-ready summaries. In 819 macaque brain samples, leave-one-sample-out Network Top1 and Top3 accuracy were 55.8% and 88.0%; leave-one-monkey-out accuracy was 53.2% and 86.7%. External analyses retained partial coarse anatomical information in normal human brain and high broad-anatomy Top3 candidate coverage in paired glioma transcriptome-MRI data. Biofluid analyses are explicitly reported as transfer stress tests when patient-level anatomical truth is unavailable. cfRNA-BrainTrace therefore provides a reproducible implementation for coarse candidate ranking and for auditing the resolution limits of atlas-based brain RNA source tracing.

### Availability and Implementation

Python 3.11+; command-line and Streamlit interfaces. Source code, documentation, tests, example data and a versioned release will be available at **[PUBLIC GITHUB OR ARCHIVAL URL REQUIRED BEFORE SUBMISSION]**. A live demonstration will be available at **[STABLE APPLICATION URL REQUIRED BEFORE SUBMISSION]**. The software will be released under **[OSI-APPROVED LICENCE REQUIRED BEFORE SUBMISSION]**.

### Contact

**[CORRESPONDING AUTHOR EMAIL REQUIRED BEFORE SUBMISSION]**

### Supplementary information

Supplementary data are available online.

## 1 Introduction

Cell-free RNA can contain recoverable tissue-of-origin information, but anatomical inference within the brain is difficult because regions share transcriptional programs and biofluid RNA is diluted, sparse and mixed with non-brain sources (Vorperian *et al.*, 2022). Existing expression-atlas workflows often return a single label without exposing whether the data support that resolution. This is problematic for translational studies, where a broad candidate may be reproducible even when exact-region localization is not.

cfRNA-BrainTrace implements a hierarchical alternative. It uses the Bo2023 macaque brain transcriptomic atlas (Bo *et al.*, 2023) to rank source candidates at four levels: lobe, broad anatomy, a 10-class functional-anatomical Network and exact region. The application treats Network and coarse Top3 candidates as the principal outputs, while exposing low-margin, low-coverage and out-of-scope conditions. It is designed for researchers who need a reproducible scoring interface, standardized diagnostics and explicit claim boundaries across tissue RNA-seq and prospective liquid-biopsy studies.

## 2 Software description

### 2.1 Inputs and preprocessing

The application accepts a two-column gene-expression table containing gene symbols and TPM-like abundance values; optional sample, subject, diagnosis and anatomical metadata are retained in exported reports. Gene identifiers are normalized and intersected with the selected model markers. External TPM-like inputs are transformed with `log1p` before scoring. This operation is a cross-scale correlation procedure, not a reconstruction of DESeq2 variance-stabilized values. The software reports matched-marker counts and non-zero coverage so that sparse inputs can be distinguished from well-covered tissue profiles.

### 2.2 Hierarchical inference

The locked model was trained from 819 samples spanning 110 macaque brain regions and nine individuals. Within each validation fold, 200 discriminative genes were selected and class centroids were scored by Pearson correlation. A pairwise rescue model re-evaluates the three leading Networks when the Top1-Top2 correlation margin is at most 0.002. The rescue threshold is stored with the model and applied identically by the command-line and web interfaces.

Predictions are propagated to lobe and broad-anatomy groupings through a versioned anatomical dictionary. Exact-region output is exploratory and is constrained by the candidate Network rather than presented as an independent clinical endpoint. The current 10-Network reference excludes cerebellum; posterior-fossa or cerebellar samples are therefore marked out of scope rather than mapped to the nearest available class.

### 2.3 Interfaces and outputs

The Streamlit interface supports file upload, metadata review, ranked Network and hierarchical candidate display, downloadable tables and model warnings. The command-line interface supports scripted scoring, benchmark export and reproducible result bundles. Each sample report contains Top1 and Top3 candidates, correlation scores, softmax-normalized display probabilities, Top1-Top2 margin, normalized entropy, marker coverage, rescue status and resolution-specific labels. Batch exports preserve model version, input settings and source-data tables.

The application separates prediction from validation. Accuracy is calculated only when an independently defined anatomical truth label is supplied. Samples without truth are summarized by candidate distribution, confidence, entropy, stability and coverage, preventing unlabeled biofluid cohorts from being described as localization validation.

## 3 Implementation

cfRNA-BrainTrace is implemented in Python using NumPy, pandas and SciPy for scoring, scikit-learn utilities for evaluation, Plotly and Matplotlib for visualization, and Streamlit for the web interface. Models are distributed as versioned NumPy and JSON artifacts containing marker order, class centroids, anatomical mappings and rescue parameters. The production Network route is implemented in `core/network_tracing.py`; the CLI and web application call the same scoring functions to avoid interface-specific divergence.

The repository includes unit tests for model loading, marker alignment, finite-value handling, probability normalization and deterministic ranking. Benchmark scripts reproduce the source-domain validation and export figures with their underlying CSV data. Optional exploratory domain-adaptation code is retained separately and is not invoked by the default scorer because its performance was endpoint-dependent and partly reliant on target-cohort information (Supplementary Methods).

## 4 Validation

Strict leave-one-sample-out validation yielded Network Top1 accuracy of 55.8% and Top3 accuracy of 88.0%. Donor-isolated leave-one-monkey-out validation yielded 53.2% and 86.7%, respectively, indicating that the main ranking performance was not explained solely by repeated sampling from the same individuals (Fig. 1B).

Cross-species evaluation used harmonized labels from the Allen Human Brain Atlas (Hawrylycz *et al.*, 2012). Among 233 samples with supported coarse labels, Network Top1 and Top3 accuracy were 32.6% and 55.4%, and lobe Top1 accuracy was 44.2%. Fine-level transfer was weaker, supporting the software's hierarchical reporting policy.

In 65 patients with paired TCGA-LGG RNA-seq and BraTS tumor segmentation (Bakas *et al.*, 2017), broad-anatomy Top3 strict and tolerant coverage were 75.4% and 83.1%; Network Top3 was 21.9% and 35.9% among 64 in-scope patients (Fig. 1C). These results support coarse candidate ranking in tumor tissue but not single-region localization.

Three public liquid-biopsy cohorts comprising serum extracellular-vesicle, tumor-associated extracellular-vesicle and cerebrospinal-fluid RNA were processed as external transfer stress tests. Predictions were generally low-margin, high-entropy or preprocessing-sensitive, despite 67/67 numerical implementation checks passing in an independent TPM-versus-count-derived CPM audit. Because these cohorts lacked patient-level imaging truth, no localization accuracy was calculated. Full validation designs, confidence intervals, disease-domain analyses and negative adaptation results are provided in the Supplementary Material.

## 5 Conclusion

cfRNA-BrainTrace packages a validated atlas-correlation workflow into reusable command-line and web interfaces and makes anatomical resolution an explicit part of prediction. Its principal use is coarse brain-origin candidate ranking with transparent confidence and coverage diagnostics. It does not convert unlabeled biofluid predictions into clinical localization claims. Future releases require human biofluid background models, calibrated abstention rules and independent cohorts with patient-level anatomical truth.

## Funding

**[FUNDING INFORMATION TO BE INSERTED BEFORE SUBMISSION]**

## Conflict of Interest

None declared.

## Data availability

The Bo2023 macaque transcriptomic atlas, Allen Human Brain Atlas, TCGA, BraTS, Ivy Glioblastoma Atlas Project and GEO datasets are available from their original repositories under the access conditions described in the Supplementary Material. Processed evaluation tables required to reproduce the reported results will be included in the versioned software release.

## References

Bakas,S. *et al.* (2017) Advancing The Cancer Genome Atlas glioma MRI collections with expert segmentation labels and radiomic features. *Sci. Data*, **4**, 170117.

Bo,T. *et al.* (2023) Brain-wide and cell-specific transcriptomic insights into MRI-derived cortical morphology in macaque monkeys. *Nat. Commun.*, **14**, 1283.

Hawrylycz,M.J. *et al.* (2012) An anatomically comprehensive atlas of the adult human brain transcriptome. *Nature*, **489**, 391-399.

Vorperian,S.K. *et al.* (2022) Cell types of origin of the cell-free transcriptome. *Nat. Biotechnol.*, **40**, 855-861.

## Figure legend

**Fig. 1. cfRNA-BrainTrace workflow and validation.** (A) Software inputs, locked scoring route and hierarchical outputs. The same scoring core is used by the command-line and Streamlit interfaces. (B) Network Top1 and Top3 accuracy in strict leave-one-sample-out (LOSO) and leave-one-monkey-out (LOMO) validation. (C) External coarse-resolution performance in normal human brain and paired TCGA-LGG/BraTS data. AHBA reports harmonized Network accuracy; paired glioma reports strict and tolerant Top3 coverage. (D) Interpretation policy. Tissue data with independent truth can support accuracy estimates, whereas unlabeled EV-RNA and CSF-RNA inputs are restricted to transfer, confidence and stability diagnostics.

## Table

**Table 1. Principal software functions and outputs.** See `tables_application_note_20260614/Table1_cfRNA_BrainTrace_features.md`.
