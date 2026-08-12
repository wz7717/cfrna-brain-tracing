# P0-3 donor-blocked DESeq2 marker audit

This analysis provides formal differential-expression support for the locked BrainTrace Network marker panel without changing the production model.

## Primary prespecified design

- Input: the immutable Bo2023 819-sample raw feature-count matrix, summed within each observed donor-by-Network cell to avoid treating multiple tissue dissections from one animal as independent biological replication.
- Biological unit and blocking factor: donor-by-Network pseudobulk, with macaque donor (`donor_id`, 9 levels) as a block.
- Anatomical factor: canonical BrainTrace Network (`network`, 10 levels); the two documented region-level corrections are applied without modifying the source workbook.
- Primary model: DESeq2 negative-binomial likelihood-ratio test, full `~ donor_id + network` versus reduced `~ donor_id`.
- Prefilter: pseudobulk raw count at least 10 in at least 3 donor-by-Network units, declared before model fitting.
- Primary multiplicity control: Benjamini-Hochberg across all tested genes for the global Network-effect LRT.
- Effect-direction follow-up: DESeq2 size-factor-normalized pseudobulk log2 counts, summarized as the highest-Network mean minus the mean of the other Networks. This is descriptive and does not create another test family.
- Locked-panel audit: map the frozen 200 Network genes to the global LRT and report the tested/mapped fraction, FDR-supported fraction and rank concordance with the historical Fisher-like score.

The pseudobulk LRT is the primary inferential audit. The 819-tissue donor-fixed-effect LRT is retained only as a sensitivity analysis because it can still overstate precision when multiple anatomical dissections from one animal are treated as residual observation units. Neither analysis replaces the fold-local predictive marker-selection rule, establishes cross-species or cfRNA validity, or provides new held-out localization accuracy.

## Reproduction

The source counts and metadata are not redistributed. Run `prepare_metadata.py` with local authorized Bo2023 paths, then run `run_deseq2_pseudobulk_marker_audit.R` for the primary audit. `run_deseq2_marker_audit.R` reproduces the tissue-level sensitivity analysis. The generated `input_manifest.json` records source hashes, and output hashes are frozen in `SHA256SUMS.txt` after a successful run.
