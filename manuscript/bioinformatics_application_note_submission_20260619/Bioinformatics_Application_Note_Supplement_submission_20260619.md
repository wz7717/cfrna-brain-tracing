# Supplementary Material for cfRNA-BrainTrace

## Supplementary Methods

### S1 Formal validation route

All validation analyses followed the same three-tier architecture used by the software. A query profile was first represented in logCPM- or logTPM-compatible expression space and then projected into Bo2023-like VSD space for Network scoring only. The top three Networks formed the candidate beam. Resolution-group and exact-region candidates were then reranked within that beam using local logCPM-compatible expression. Fold-local Network model construction differed between the LOSO and LOMO validation designs as documented in S7, but projected VSD remained restricted to Network beam generation in both. This design separates broad candidate generation from fine regional interpretation.

### S2 Internal validation design

Internal validation used the Bo2023 macaque brain reference. Two settings were evaluated. Leave-one-sample-out validation tested whether the route could recover the held-out sample label when other samples from the reference remained available. Leave-one-monkey-out validation held out all samples from one animal at a time and therefore tested donor-level generalization. Network metrics included all 819 samples because every held-out sample retained a supported Network label. Resolution-group and exact-region metrics required the truth region to remain represented in the corresponding training fold. Five LOSO samples and seven LOMO samples did not meet that region-reference requirement and were excluded only from region-level denominators. Top1, Top3 and median true-rank were used to distinguish exact calls from candidate-list recovery.

### S3 External validation design

External validation was limited by the label resolution available in each dataset. AHBA human brain RNA-seq was used for mapped-label validation because its anatomical labels could be harmonized to Network, resolution-group and a subset of exact-region labels. TCGA/BraTS glioma tissue RNA-seq with MRI-derived labels was used only for coarse anatomical consistency because its truth labels are human imaging labels, not Bo2023 macaque exact-region identifiers. GSE189919 was used to test whether an external matrix could be projected into the model gene space; it was not used for accuracy estimation because patient-level anatomical truth was unavailable.

### S4 Label harmonization and allowed conclusions

Label harmonization was performed only to the level supported by each dataset. AHBA anatomical labels were mapped to the macaque-derived Network, resolution-group and exact-region hierarchy where a supported mapping existed, and results are interpreted as mapped-label transfer rather than direct anatomical equivalence. TCGA/BraTS labels were derived from human MRI/tumour context and therefore support only coarse tumour-tissue anatomical consistency. Biofluid datasets without patient-level anatomical truth were not used for localization accuracy and are reported only as projection-feasibility or transfer stress tests.

### S5 Model artifacts and projected-VSD projector

The repository contains lightweight model artifacts in `data/models/`, including Network model files, region-resolution dictionaries, route metadata and the reference-fitted projector `bo2023_reference_projector_linear_full.npz`. The projector stores gene-wise slope, intercept and clipping parameters fitted from reference data. At inference time, a single query profile can be mapped into Bo2023-like projected-VSD space using these fixed parameters without target-cohort labels or target-cohort distribution fitting. Projected VSD is used only for Network Top3 beam generation; resolution-group and exploratory exact-region reranking use logCPM-compatible local expression within that beam.

### S6 Example input/output

Synthetic public example files are provided in `submission_ready_assets/example_io/` as part of the v0.1.6 public submission release. They document accepted input columns, Network ranked candidates, three-tier JSON output, resolution-group ranked candidates, exact-region ranked candidates and the local generation script. The example files are format and reproducibility aids only; they are not biological validation samples.

### S6a Software availability

The manuscript-associated software release is v0.1.6. The public GitHub repository is available at `https://github.com/wz7717/cfrna-brain-tracing`, and the archived release is cited in the manuscript using the Zenodo version DOI `https://doi.org/10.5281/zenodo.20780280`; the project concept DOI is `https://doi.org/10.5281/zenodo.20773674`.

### S7 Model development and locked evaluation timeline

The three-tier architecture was fixed before generating the final submission validation tables. External AHBA, TCGA/BraTS and biofluid analyses were not used to select the final route; they were used only for mapped-label transfer evaluation, coarse tumour-tissue consistency assessment and biofluid transfer stress testing. The LOSO implementation used locked Network genes and correlation-ranked projected-VSD Network Top3 candidates. The formal LOMO implementation rebuilt discriminative Network genes and a fold-local pairwise Top1 rescue within each training fold; this rescue could reorder Top1 but did not change the original Network Top3 candidate set. Both implementations used projected VSD only for the Network beam and logCPM-compatible local evidence for downstream reranking. The independent Network-only LOSO and LOMO analyses are reported separately as route-selection evidence.

### S8 Figures and tables

Figure source data and synthetic example input/output files are archived in the GitHub/Zenodo release rather than treated as separate journal supplementary files. Supplementary Tables S1-S6 are embedded below so that the formal supplementary PDF can be submitted as a single file containing Supplementary Methods, Supplementary Results and Tables S1-S6. Together, these materials document both how the route was tested and what result each validation setting supports.

## Supplementary Results

### S1 Internal route selection

The first internal analysis evaluated Network-level candidate generation. Projected VSD achieved Network Top1/Top3 of 58.00%/91.58% in LOSO and 53.72%/91.33% in LOMO, exceeding logCPM baseline and native VSD at Network Top3. Direct exact-region scoring was lower and less stable, especially in LOMO. These results motivated the use of projected VSD for broad Network beam generation and logCPM-compatible expression for downstream local reranking.

### S2 Formal internal three-tier validation

In the complete LOSO validation, the formal route achieved Network Top1/Top3 of 58.24%/92.19% across all 819 samples. Resolution-group and exact-region Top3 were 72.36% and 45.33% among 814 reference-supported samples. In complete LOMO validation, Network Top3 was 91.21% across all 819 samples, while resolution-group and exact-region Top3 were 69.09% and 42.36% among 812 reference-supported samples. The 92.19% LOSO Network Top3 value uses all 819 Network-evaluable samples as the denominator. Region-level LOSO metrics use 814 reference-supported samples because five samples lacked a truth-region reference after fold construction. Region-level LOMO metrics use 812 reference-supported samples because seven samples lacked a truth-region reference after fold construction. The earlier LOSO Network value of 92.38% was conditional on the 814 region-evaluable samples and is retained only as a legacy denominator inconsistency, not as the submission result. Median true-rank increased from Network to exact-region levels, consistent with decreasing anatomical certainty at finer resolution. Resolution group is therefore the preferred region-level endpoint, while exact-region output is retained as a candidate ranking.

### S3 External validation results

In AHBA, the hybrid route achieved Network Top1/Top3 of 74.68%/94.42%, resolution-group Top1/Top3 of 36.26%/67.03% and exact-region Top1/Top3 of 24.18%/42.86%. Hybrid exact Top3 exceeded logCPM baseline (30.77%) and projected-VSD-only scoring (29.67%). In TCGA/BraTS, hybrid Network Top3 was 40.00% and broad-anatomy Top3 was 64.62%, supporting only coarse anatomical consistency. GSE189919 overlapped 15,622/21,668 projector genes, corresponding to 72.10% gene-space coverage, supporting projection feasibility rather than source-localization accuracy.

## Supplementary Tables

### Table S1. Internal validation design

| Validation | Data | Held-out unit | Route tested | Reported endpoints | Denominator policy |
|---|---|---|---|---|---|
| Network LOSO | Bo2023 macaque brain RNA-seq | Single sample | Projected VSD vs logCPM/native VSD Network scoring | Network Top1, Network Top3, median true-rank | All 819 samples |
| Network LOMO | Bo2023 macaque brain RNA-seq | One monkey | Projected VSD vs logCPM/native VSD Network scoring | Network Top1, Network Top3, median true-rank | All 819 samples |
| Formal three-tier LOSO | Bo2023 macaque brain RNA-seq | Single sample | Projected-VSD Network beam plus logCPM group/exact rerank | Network, resolution-group and exact-region Top1/Top3 | Network n=819; group/exact n=814 |
| Formal three-tier LOMO | Bo2023 macaque brain RNA-seq | One monkey | Fold-local Network beam plus logCPM group/exact rerank | Network, resolution-group and exact-region Top1/Top3 | Network n=819; group/exact n=812 |

### Table S2. Internal validation results

| Dataset | Route | Endpoint | Evaluated n | Top1 | Top3 | Median true-rank | Interpretation |
|---|---|---|---:|---|---|---|---|
| Bo2023 LOSO | Projected VSD Network | Network | 819 | 58.00% | 91.58% | 1.0 | Supports projected-VSD Network beam |
| Bo2023 LOMO | Projected VSD Network | Network | 819 | 53.72% | 91.33% | 1.0 | Supports donor-level Network beam |
| Bo2023 LOSO | Formal hybrid | Network | 819 | 58.24% | 92.19% | 1.0 | Primary endpoint across all Network-evaluable samples |
| Bo2023 LOSO | Formal hybrid | Resolution group | 814 | 44.47% | 72.36% | 2.0 | Main region-level endpoint among reference-supported samples |
| Bo2023 LOSO | Formal hybrid | Exact region | 814 | 22.48% | 45.33% | 4.0 | Exploratory ranking among reference-supported samples |
| Bo2023 LOMO | Formal hybrid | Network | 819 | 57.75% | 91.21% | 1.0 | Cross-monkey support |
| Bo2023 LOMO | Formal hybrid | Resolution group | 812 | 41.38% | 69.09% | 2.0 | Main region-level endpoint among reference-supported samples |
| Bo2023 LOMO | Formal hybrid | Exact region | 812 | 22.17% | 42.36% | 5.0 | Exploratory ranking among reference-supported samples |

### Table S3. External validation design

| Dataset | Sample type | n | Truth label type | Allowed conclusion |
|---|---|---|---|---|
| AHBA | Human normal brain RNA-seq | 242 total; 233 supported; 91 exact-evaluable | Mapped anatomical labels | Cross-species mapped-label validation |
| TCGA/BraTS | Glioma tissue RNA-seq with MRI-derived labels | 65 patients | Human imaging labels | Coarse anatomical consistency only |
| GSE189919 | External count matrix | 51 samples | No patient-level anatomical truth | Projection feasibility only |

### Table S4. External validation results

| Dataset | Route | Endpoint | Top1 | Top3 | Conclusion |
|---|---|---|---|---|---|
| AHBA | Hybrid | Network | 74.68% | 94.42% | Strong mapped-label Network support |
| AHBA | Hybrid | Resolution group | 36.26% | 67.03% | Moderate group-level support |
| AHBA | Hybrid | Exact region | 24.18% | 42.86% | Only exact-evaluable mapped labels |
| AHBA | logCPM baseline | Exact region | 17.58% | 30.77% | Below hybrid exact Top3 |
| AHBA | Projected VSD only | Exact region | 10.99% | 29.67% | Below hybrid exact Top3 |
| TCGA/BraTS | Hybrid | Network | 15.38% | 40.00% | Coarse consistency only |
| TCGA/BraTS | Hybrid | Broad anatomy | 13.85% | 64.62% | Coarse consistency only |
| GSE189919 | Projection feasibility | Projector gene overlap | 15622/21668 | 72.10% | No accuracy claim |

### Table S5. Figure and table index

| Item | Content | Use in manuscript |
|---|---|---|
| Figure 1 | Two-panel workflow and validation summary; Panel B numeric source in `Figure1_validation_summary.csv` | Main Application Note figure |
| Table S1 | Internal validation design | Documents LOSO/LOMO operations |
| Table S2 | Internal validation results | Reports Network, group and exact-region metrics |
| Table S3 | External validation design | Documents label support and allowed conclusions |
| Table S4 | External validation results | Reports AHBA, TCGA/BraTS and GSE189919 outcomes |
| Table S5 | Figure/table index | Indexes manuscript figure and supplementary tables |
| Table S6 | Claim boundaries | Defines unsupported interpretations and permitted claims |

### Table S6. Claim boundaries

| Avoid | Use |
|---|---|
| Projected VSD creates a new Bo2023 atlas. | Only query profiles are projected for Network beam generation. |
| Projected VSD is best for exact-region inference. | Projected VSD supports Network beam; logCPM supports downstream reranking. |
| TCGA/BraTS validates Bo2023 exact regions. | TCGA/BraTS supports coarse anatomical consistency. |
| GSE189919 validates accuracy. | GSE189919 verifies projection feasibility only. |
| Exact Top1 is deterministic localization. | Exact outputs are exploratory local candidate rankings. |
