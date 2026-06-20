# Supplementary Material for cfRNA-BrainTrace

## Supplementary Methods

### S1 Formal validation route

All validation analyses followed the same formal route used by the software. A query profile was first represented in logCPM- or logTPM-compatible expression space and then projected into Bo2023-like VSD space for Network scoring only. The top three Networks formed the candidate beam. Resolution-group and exact-region candidates were then reranked within that beam using local logCPM-compatible expression. This design separates broad candidate generation from fine regional interpretation.

### S2 Internal validation design

Internal validation used the Bo2023 macaque brain reference. Two settings were evaluated. Leave-one-sample-out validation tested whether the route could recover the held-out sample label when other samples from the reference remained available. Leave-one-monkey-out validation held out all samples from one animal at a time and therefore tested donor-level generalization. Metrics were reported separately for Network, resolution-group and exact-region levels. Top1, Top3 and median true-rank were used to distinguish exact calls from candidate-list recovery.

### S3 External validation design

External validation was limited by the label resolution available in each dataset. AHBA human brain RNA-seq was used for mapped-label validation because its anatomical labels could be harmonized to Network, resolution-group and a subset of exact-region labels. TCGA/BraTS glioma tissue RNA-seq with MRI-derived labels was used only for coarse anatomical consistency because its truth labels are human imaging labels, not Bo2023 macaque exact-region identifiers. GSE189919 was used to test whether an external matrix could be projected into the model gene space; it was not used for accuracy estimation because patient-level anatomical truth was unavailable.

### S4 Label harmonization and allowed conclusions

Label harmonization was performed only to the level supported by each dataset. AHBA anatomical labels were mapped to the macaque-derived Network, resolution-group and exact-region hierarchy where a supported mapping existed, and results are interpreted as mapped-label transfer rather than direct anatomical equivalence. TCGA/BraTS labels were derived from human MRI/tumour context and therefore support only coarse tumour-tissue anatomical consistency. Biofluid datasets without patient-level anatomical truth were not used for localization accuracy and are reported only as projection-feasibility or transfer stress tests.

### S5 Model artifacts and projected-VSD projector

The repository contains lightweight model artifacts in `data/models/`, including Network model files, region-resolution dictionaries, route metadata and the reference-fitted projector `bo2023_reference_projector_linear_full.npz`. The projector stores gene-wise slope, intercept and clipping parameters fitted from reference data. At inference time, a single query profile can be mapped into Bo2023-like projected-VSD space using these fixed parameters without target-cohort labels or target-cohort distribution fitting. Projected VSD is used only for Network Top3 beam generation; resolution-group and exploratory exact-region reranking use logCPM-compatible local expression within that beam.

### S6 Example input/output

Synthetic public example files are provided in `submission_ready_assets/example_io/`. They document accepted input columns, Network ranked candidates, three-tier JSON output, resolution-group ranked candidates, exact-region ranked candidates and the local generation script. The example files are format and reproducibility aids only; they are not biological validation samples.

### S7 Model development and locked evaluation timeline

The production route was fixed before generating the final submission validation tables. External AHBA, TCGA/BraTS and biofluid analyses were not used to select the final route; they were used only for mapped-label transfer evaluation, coarse tumour-tissue consistency assessment and biofluid transfer stress testing. Development comparisons that used the same internal validation framework are reported as development evidence rather than independent confirmation.

### S8 Figures and tables

Supplementary Figure S1 provides the main route and validation summary artwork. Supplementary Tables S1-S6 provide the internal validation design, internal validation results, external validation design, external validation results, figure/table index and claim-boundary summary. Together, these materials document both how the route was tested and what result each validation setting supports.

## Supplementary Results

### S1 Internal route selection

The first internal analysis evaluated Network-level candidate generation. Projected VSD achieved Network Top1/Top3 of 58.00%/91.58% in LOSO and 53.72%/91.33% in LOMO, exceeding logCPM baseline and native VSD at Network Top3. Direct exact-region scoring was lower and less stable, especially in LOMO. These results motivated the use of projected VSD for broad Network beam generation and logCPM-compatible expression for downstream local reranking.

### S2 Formal internal three-tier validation

In the complete LOSO validation, the formal route achieved Network, resolution-group and exact-region Top3 values of 92.38%, 72.36% and 45.33%. In complete LOMO validation, the corresponding Top3 values were 91.21%, 69.09% and 42.36%. Median true-rank increased from Network to exact-region levels, consistent with decreasing anatomical certainty at finer resolution. Resolution group is therefore the preferred region-level endpoint, while exact-region output is retained as a candidate ranking.

### S3 External validation results

In AHBA, the hybrid route achieved Network Top1/Top3 of 74.68%/94.42%, resolution-group Top1/Top3 of 36.26%/67.03% and exact-region Top1/Top3 of 24.18%/42.86%. Hybrid exact Top3 exceeded logCPM baseline (30.77%) and projected-VSD-only scoring (29.67%). In TCGA/BraTS, hybrid Network Top3 was 40.00% and broad-anatomy Top3 was 64.62%, supporting only coarse anatomical consistency. GSE189919 overlapped 15,622/21,668 projector genes, corresponding to 72.10% gene-space coverage, supporting projection feasibility rather than source-localization accuracy.
