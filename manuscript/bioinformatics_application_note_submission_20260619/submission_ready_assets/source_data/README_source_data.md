# Source Data Working Package

This directory is prepared for public archival with the manuscript-associated
v0.1.6 release.
It contains local copies of the current supplementary tables and Figure 1 artwork.

## Included Files

- `TableS1_internal_validation_design.csv`
- `TableS2_internal_validation_results.csv`
- `TableS3_external_validation_design.csv`
- `TableS4_external_validation_results.csv`
- `TableS5_figure_table_index.csv`
- `TableS6_claim_boundaries.csv`
- `Figure1_validation_summary.csv`
- `Figure1_cfRNA_BrainTrace_final.svg`
- `Figure1_cfRNA_BrainTrace_final.pdf`
- `Figure1_cfRNA_BrainTrace_final.png`
- `Figure1_cfRNA_BrainTrace_Bioinformatics_lowres.png`
- `Figure1_cfRNA_BrainTrace_Bioinformatics_highres_178mm.tif`

## Status

Archived with the manuscript-associated v0.1.6 GitHub/Zenodo release. The manuscript
uses the Zenodo version DOI `https://doi.org/10.5281/zenodo.20780280`. The project
concept DOI is `https://doi.org/10.5281/zenodo.20773674`.

## Notes

The tables are derived summary/evaluation tables intended for manuscript support.
`Figure1_validation_summary.csv` is the authoritative numeric source for the
Figure 1 Panel B validation bars. The final Figure 1 artwork is a two-panel
workflow and validation summary, with vector source files preserved as SVG/PDF.
Network LOSO uses all 819 samples; region-level
LOSO uses 814 reference-supported samples; region-level LOMO uses 812
reference-supported samples. The 92.19% LOSO Network Top3 value uses all 819
Network-evaluable samples as the denominator. Five LOSO samples and seven LOMO
samples lacked a truth-region reference after fold construction and are excluded
only from region-level evaluation.
They do not include patient-level MRI, clinical data, raw TCGA/GEO/AHBA matrices, or raw Bo2023 source matrices.
