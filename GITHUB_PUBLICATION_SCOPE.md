# GitHub Publication Scope

## Included in the first repository upload

- Python application, CLI and core inference code.
- Lightweight locked model artifacts and anatomical dictionaries.
- Tests and environment/configuration files.
- Validation, figure and manuscript build scripts.
- Bioinformatics Application Note drafts, figures, tables and supplementary PDF.
- Documentation required to understand installation, inputs and claim boundaries.

## Excluded

- `cfrna_source_tracing.db`.
- Raw or downloaded data under `data/`, except `data/models/`.
- `bo2023 data/`, TCGA, BraTS, MRI, GEO, AHBA and other source matrices.
- `results/`, local reports and generated analysis caches.
- `vendor/` and `tools/`, including third-party executables.
- Streamlit secrets, credentials and environment files.
- Temporary render directories and local package installations.

## Before making the repository public

1. Select and add an OSI-approved software license.
2. Confirm that every model and derived annotation can be redistributed.
3. Replace manuscript URL, author, contact and funding placeholders.
4. Create a tagged release and archive it with a DOI.
5. Verify installation and tests from a clean checkout without private data.
