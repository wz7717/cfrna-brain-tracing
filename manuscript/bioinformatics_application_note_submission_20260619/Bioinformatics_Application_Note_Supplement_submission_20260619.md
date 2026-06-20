# Supplementary Material for cfRNA-BrainTrace

## Supplementary Methods

### S1 Formal validation route

All validation analyses followed the same formal route used by the software. A query profile was first represented in logCPM- or logTPM-compatible expression space and then projected into Bo2023-like VSD space for Network scoring only. The top three Networks formed the candidate beam. Resolution-group and exact-region candidates were then reranked within that beam using local logCPM-compatible expression. This design separates broad candidate generation from fine regional interpretation.

### S2 Internal validation design

Internal validation used the Bo2023 macaque brain reference. Two settings were evaluated. Leave-one-sample-out validation tested whether the route could recover the held-out sample label when other samples from the reference remained available. Leave-one-monkey-out validation held out all samples from one animal at a time and therefore tested donor-level generalization. Metrics were reported separately for Network, resolution-group and exact-region levels. Top1, Top3 and median true-rank were used to distinguish exact calls from candidate-list recovery.

### S3 External validation design

External validation was limited by the label resolution available in each dataset. AHBA human brain RNA-seq was used for mapped-label validation because its anatomical labels could be harmonized to Network, resolution-group and a subset of exact-region labels. TCGA/BraTS glioma tissue RNA-seq with MRI-derived labels was used only for coarse anatomical consistency because its truth labels are human imaging labels, not Bo2023 macaque exact-region identifiers. GSE189919 was used to test whether an external matrix could be projected into the model gene space; it was not used for accuracy estimation because patient-level anatomical truth was unavailable.

### S4 Figures and tables

Supplementary Figure S1 provides the main route and validation summary artwork. Supplementary Tables S1-S5 provide the internal validation design, internal validation results, external validation design, external validation results and figure/table index. Together, these materials document both how the route was tested and what result each validation setting supports.

## Supplementary Results

### S1 Internal route selection

The first internal analysis evaluated Network-level candidate generation. Projected VSD achieved Network Top1/Top3 of 58.00%/91.58% in LOSO and 53.72%/91.33% in LOMO, exceeding logCPM baseline and native VSD at Network Top3. Direct exact-region scoring was lower and less stable, especially in LOMO. These results motivated the use of projected VSD for broad Network beam generation and logCPM-compatible expression for downstream local reranking.

### S2 Formal internal three-tier validation

In the complete LOSO validation, the formal route achieved Network, resolution-group and exact-region Top3 values of 92.38%, 72.36% and 45.33%. In complete LOMO validation, the corresponding Top3 values were 91.21%, 69.09% and 42.36%. Median true-rank increased from Network to exact-region levels, consistent with decreasing anatomical certainty at finer resolution. Resolution group is therefore the preferred region-level endpoint, while exact-region output is retained as a candidate ranking.

### S3 External validation results

In AHBA, the hybrid route achieved Network Top1/Top3 of 74.68%/94.42%, resolution-group Top1/Top3 of 36.26%/67.03% and exact-region Top1/Top3 of 24.18%/42.86%. Hybrid exact Top3 exceeded logCPM baseline (30.77%) and projected-VSD-only scoring (29.67%). In TCGA/BraTS, hybrid Network Top3 was 40.00% and broad-anatomy Top3 was 64.62%, supporting only coarse anatomical consistency. GSE189919 overlapped 15,622/21,668 projector genes, corresponding to 72.10% gene-space coverage, supporting projection feasibility rather than source-localization accuracy.
