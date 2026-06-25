# cfRNA-BrainTrace v0.1.6

This is the cleaned public submission release for the Bioinformatics Application Note describing cfRNA-BrainTrace.

## Update

- Merged the public submission cleanup into the default `main` branch.
- Removed manuscript-irrelevant historical drafts, internal handoff materials, private atlas/API deployment files, obsolete figure/table build outputs and development-only audit/download scripts from the public repository.
- Retained the final English manuscript, supplementary files, Tables S1-S6, Figure 1, figure source data, synthetic example input/output, app code, core models, tests and current validation/benchmark materials.
- Synchronized README, manuscript and supplementary availability statements to the manuscript-associated v0.1.6 release.
- Current submission-route validation numbers include the denominator correction documented below.

## Validation denominator correction

- Formal LOSO Network metrics now include all 819 Network-evaluable samples:
  Top1 58.24% and Top3 92.19%.
- Resolution-group and exact-region LOSO metrics remain restricted to 814
  reference-supported samples; five held-out truth regions are absent from
  their training folds.
- Formal LOMO Network metrics use all 819 samples, while region-level metrics
  use 812 reference-supported samples.
- The former formal LOSO Network Top3 value of 92.38% is a legacy conditional
  result on the 814 region-evaluable samples and is not the submission value.

## Archive identifiers

- Version DOI: https://doi.org/10.5281/zenodo.20780280
- Concept DOI: https://doi.org/10.5281/zenodo.20773674
