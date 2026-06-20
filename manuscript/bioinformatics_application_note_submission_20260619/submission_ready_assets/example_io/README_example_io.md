# Public Example Input/Output Package

This directory contains a synthetic, redistributable example input and corresponding output files generated with the project core scoring functions.

## Included

- `example_expression_input_template.csv`: small format-only input template.
- `example_expression_input_full_200genes.csv`: synthetic TPM-like expression values for all 200 Network model marker genes.
- `example_network_output.json`: full JSON output from Network scoring.
- `example_network_ranked_candidates.csv`: ranked Network candidates.
- `example_three_tier_output.json`: full JSON output from the three-tier route.
- `example_resolution_group_ranked_candidates.csv`: resolution-group candidate table derived from the three-tier output.
- `example_exact_region_ranked_candidates.csv`: exact-region candidate table from the three-tier output.
- `generate_example_output.py`: local reproduction script for the example outputs.

## Important Boundary

The example is synthetic and is intended to demonstrate file format and output structure only.
It is not a biological validation sample and should not be described as evidence of tracing performance.

## Generation Status

- Dependencies were installed in the active Python environment.
- `cli.py --help` runs successfully after dependency installation.
- Network scoring on `example_expression_input_full_200genes.csv` produced high traceability with 200/200 model-gene overlap.
- Three-tier scoring produced high traceability and 30 exact-region candidates.
- No patient-level, clinical, MRI, GEO, TCGA, AHBA or raw Bo2023 data are included.

## Reproducibility Note

The generated outputs use `core.network_tracing.trace_network_expression` and `core.bo2023_region_tracing.trace_bo2023_secondary_regions` with the local model artifacts and reference resources available in the repository/workspace.

From the repository root, regenerate the example outputs with:

```bash
python manuscript/bioinformatics_application_note_submission_20260619/submission_ready_assets/example_io/generate_example_output.py
```
