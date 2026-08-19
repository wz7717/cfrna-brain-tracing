# Retired Huang 2025 pseudo-pairing material

## Status

**NOT USED FOR MANUSCRIPT**

**INVALID FOR PATIENT-PAIRED INFERENCE**

**RETAINED ONLY FOR HISTORICAL TRACEABILITY**

This material is not used for manuscript claims, supplementary claims,
canonical outputs, or patient-paired inference.

Earlier local Huang 2025 analyses treated matching terminal numbers in labels
such as `GLI_CSF16` and `GLI_plasma16` as a patient-level correspondence. The
public expression matrix does not supply evidence that this label convention
identifies the same patient across fluids. Any result produced by that
assumption, including CSF-plasma concordance, within-group permutation values,
or synthetic mixture/admixture analyses, is invalid for patient-paired
interpretation.

The retired implementation and any associated local artifacts are retained
only in version-control history or separately quarantined forensic records.
The working forensic copy is isolated at
`legacy_current_script` (local historical audit workspace; not distributed) and
is not copied into `reproducibility/huang_2025/`. The canonical
replacement is the full published-matrix audit of 159 profiles analysed in
separate fluid-specific cohorts (77 CSF and 82 plasma), documented in
`../huang_2025/HUANG_2025_RESULTS.md`.

Because the public matrix does not provide a profile-to-patient map,
patient-level dependence or independence among these profiles cannot be
assessed.

This retirement does not question the source article. It corrects an
unsupported downstream assumption made in BrainTrace's prior local analysis.
