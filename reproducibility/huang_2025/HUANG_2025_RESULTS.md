# Huang 2025 cfRNA: provenance-remediated external domain audit

Source: Huang et al. (2025), DOI [10.1038/s41698-025-00909-6](https://doi.org/10.1038/s41698-025-00909-6).

## Scope and provenance

The public expression matrix was used as a computational stress-test resource. Although the source study applied its own sequencing-QC exclusions for its clinical analyses, all 159 profiles available in the published matrix were considered here for technical transfer auditing; this analysis does not attempt to reproduce the source study’s QC-filtered clinical analysis.

No patient-level CSF-plasma correspondence was assumed. CSF and plasma were analysed as separate fluid-specific profile cohorts; patient-level dependence or independence cannot be established from the public matrix. Sample-label suffixes were not interpreted as patient identifiers.

This audit supports a narrow technical-portability/domain-shift statement. It does not establish patient correspondence, synthetic mixture behavior, anatomical localization accuracy, tumour-source discrimination, or clinical validity.

## Cohort accounting

- Published-matrix audit universe: 159 profiles (77 CSF; 82 plasma).
- Traceable BrainTrace outputs: 159/159 (100.0%).
- OMPFC Top1: CSF 74/77 (96.1%); plasma 56/82 (68.3%).

## Profile-level tumour-control diagnostics

Six two-sided Mann-Whitney U tests were run across fluid and metric; the smallest Benjamini-Hochberg FDR was 0.722052. All tests are profile-level analyses with pairing unavailable; patient-level dependence cannot be verified from the public matrix.

| Fluid | Metric | tumour n | control n | raw P | BH-FDR |
|---|---|---:|---:|---:|---:|
| CSF | atlas_fit_score | 59 | 18 | 0.320726 | 0.722052 |
| CSF | network_margin | 59 | 18 | 0.481368 | 0.722052 |
| CSF | network_entropy | 59 | 18 | 0.677966 | 0.813559 |
| plasma | atlas_fit_score | 64 | 18 | 0.367151 | 0.722052 |
| plasma | network_margin | 64 | 18 | 0.968723 | 0.968723 |
| plasma | network_entropy | 64 | 18 | 0.321468 | 0.722052 |

## Exploratory marker correlations

Marker associations are descriptive, fluid-specific Spearman correlations with the OMPFC network score; they are not matched-biofluid comparisons; patient-level dependence cannot be assessed from the public matrix.

| Fluid | Marker class | Marker | n | rho | raw P | BH-FDR |
|---|---|---|---:|---:|---:|---:|
| CSF | platelet-associated | PF4 | 77 | 0.1497 | 0.193883 | 0.238625 |
| CSF | platelet-associated | PPBP | 77 | 0.0254 | 0.826658 | 0.853324 |
| CSF | platelet-associated | RGS18 | 77 | -0.0426 | 0.712833 | 0.760355 |
| CSF | platelet-associated | GP9 | 77 | -0.2755 | 0.015288 | 0.025748 |
| CSF | platelet-associated | ITGA2B | 77 | -0.0799 | 0.489963 | 0.559958 |
| CSF | platelet-associated | TUBB1 | 77 | 0.0742 | 0.521479 | 0.575425 |
| CSF | platelet-associated | SELP | 77 | -0.2428 | 0.033338 | 0.046383 |
| CSF | platelet-associated | NRGN | 77 | -0.3996 | 0.000318 | 0.000782 |
| CSF | extracellular-vesicle-associated | CD9 | 77 | -0.2099 | 0.066867 | 0.089156 |
| CSF | extracellular-vesicle-associated | CD63 | 77 | -0.4761 | 0.000012 | 0.000048 |
| CSF | extracellular-vesicle-associated | CD81 | 77 | -0.2892 | 0.010756 | 0.020246 |
| CSF | extracellular-vesicle-associated | TSG101 | 77 | -0.4119 | 0.000198 | 0.000575 |
| CSF | extracellular-vesicle-associated | PDCD6IP | 77 | -0.4019 | 0.000291 | 0.000776 |
| CSF | extracellular-vesicle-associated | SDCBP | 77 | -0.4402 | 0.000062 | 0.000198 |
| CSF | extracellular-vesicle-associated | FLOT1 | 77 | -0.3590 | 0.001347 | 0.002873 |
| CSF | extracellular-vesicle-associated | FLOT2 | 77 | -0.2516 | 0.027263 | 0.039655 |
| plasma | platelet-associated | PF4 | 82 | 0.4338 | 0.000047 | 0.000166 |
| plasma | platelet-associated | PPBP | 82 | 0.3754 | 0.000510 | 0.001165 |
| plasma | platelet-associated | RGS18 | 82 | -0.3161 | 0.003818 | 0.007636 |
| plasma | platelet-associated | GP9 | 82 | 0.0126 | 0.910802 | 0.910802 |
| plasma | platelet-associated | ITGA2B | 82 | -0.2692 | 0.014468 | 0.025721 |
| plasma | platelet-associated | TUBB1 | 82 | 0.1879 | 0.090913 | 0.116368 |
| plasma | platelet-associated | SELP | 82 | -0.7165 | 0.000000 | 0.000000 |
| plasma | platelet-associated | NRGN | 82 | 0.2529 | 0.021911 | 0.033388 |
| plasma | extracellular-vesicle-associated | CD9 | 82 | -0.6499 | 0.000000 | 0.000000 |
| plasma | extracellular-vesicle-associated | CD63 | 82 | -0.0938 | 0.401668 | 0.476052 |
| plasma | extracellular-vesicle-associated | CD81 | 82 | -0.5624 | 0.000000 | 0.000000 |
| plasma | extracellular-vesicle-associated | TSG101 | 82 | -0.7503 | 0.000000 | 0.000000 |
| plasma | extracellular-vesicle-associated | PDCD6IP | 82 | -0.7896 | 0.000000 | 0.000000 |
| plasma | extracellular-vesicle-associated | SDCBP | 82 | -0.4890 | 0.000003 | 0.000017 |
| plasma | extracellular-vesicle-associated | FLOT1 | 82 | -0.4759 | 0.000006 | 0.000028 |
| plasma | extracellular-vesicle-associated | FLOT2 | 82 | -0.2529 | 0.021901 | 0.033388 |

## Canonical outputs

The sample ledger, per-profile outputs, rankings, fluid-specific profile-cohort distributions, statistics, machine-readable summaries, and audit manifest in this directory are the canonical Huang 2025 remediation outputs. Pseudo-paired CSF-plasma and synthetic-mixture outputs are intentionally absent.
