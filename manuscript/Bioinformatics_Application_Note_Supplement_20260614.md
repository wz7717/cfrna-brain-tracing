# Supplementary Material for cfRNA-BrainTrace

## Supplementary Methods

### S1 Reference atlas and anatomical hierarchy

The Bo2023 reference contains 819 macaque brain RNA-seq samples from nine individuals, covering 110 annotated regions. Regions were mapped to lobe, broad anatomy and 10 Saleem-style Networks. The reference does not contain cerebellum; cerebellar and posterior-fossa samples are out of scope at the Network level.

### S2 Locked Network model

Within each source-domain validation fold, 200 discriminative genes were selected without using the held-out sample. Pearson correlation between the aligned sample vector and class centroids produced the primary ranking. A pairwise rescue model re-ranked the three leading classes when the Top1-Top2 margin was no greater than 0.002. The threshold was selected retrospectively and requires independent confirmation.

### S3 External input scale

External TPM-like expression is transformed with `log1p` and correlated with the Bo2023 variance-stabilized reference. Pearson correlation removes linear shifts and scaling but does not remove non-linear distribution differences, zero inflation, batch effects or sample-type effects. External analyses are therefore described as cross-scale correlation stress tests.

### S4 Anatomical truth policy

Accuracy was calculated only when labels were defined independently of expression predictions. For TCGA-LGG/BraTS, tumor masks were compared with directly overlapping valid atlas labels. Unbounded nearest-neighbour filling was excluded because it artificially expanded deep white-matter labels. Strict evaluation used the dominant direct-overlap label; tolerant evaluation accepted prespecified directly overlapping alternatives.

### S5 Exploratory adaptation

Tumor-state filtering, rank-based quantile harmonization and class-prior calibration were evaluated only as sensitivity analyses. Harmonization using target-cohort distributions was transductive; MRI-label-informed calibration was cohort-internal. The implemented rank mapping was not tie-aware and could assign equal zero-valued markers to different quantiles. These routes were excluded from the locked workflow.

## Supplementary Results

### S1 Internal validation

Strict leave-one-sample-out Network Top1 and Top3 accuracy were 55.8% and 88.0%. Leave-one-monkey-out Top1 and Top3 accuracy were 53.2% and 86.7%. Finer region-group and exact-region outputs were less stable and remain secondary or exploratory.

### S2 Normal human brain transfer

Among 233 AHBA samples with harmonized coarse labels, Network Top1 and Top3 accuracy were 32.6% and 55.4%; lobe Top1 accuracy was 44.2%. Among 91 exact-mapped samples, exact-region Top1 and Top3 accuracy were 9.9% and 29.7%. Fine localization was therefore not treated as a validated endpoint.

### S3 Paired glioma transcriptome-MRI evaluation

In 65 paired TCGA-LGG/BraTS patients, lobe Top3 strict/tolerant coverage was 84.6%/89.2% and broad-anatomy Top3 coverage was 75.4%/83.1%. Among 64 Network-evaluable patients, Network Top3 was 21.9%/35.9%. One cerebellar tumor was excluded from the Network denominator because it was outside the reference space.

### S4 Liquid-biopsy transfer stress tests

GSE228512 comprised 85 GBM and 31 healthy serum EV-RNA samples. GSE106804 comprised 13 GBM and six healthy tumor-associated EV-capture samples. GSE189919 comprised 40 medulloblastoma and 11 normal CSF-RNA samples. None contained patient-level anatomical imaging truth.

GSE228512 predictions were dominated by Subcortical outputs in both GBM and healthy samples. Tumor-associated EV capture in GSE106804 reduced single-class dominance but also reduced margins and produced only 21.1% Top1 agreement across adapted preprocessing routes. In GSE189919, all 67 numerical audit checks passed and locked-baseline Top1 agreement between TPM and count-derived CPM was 87.5% in medulloblastoma and 90.9% in controls; nevertheless, probabilities remained near uniform and normalized entropy was nearly maximal. These cohorts establish technical executability and domain mismatch, not localization accuracy.

### S5 Adaptation sensitivity

Unsupervised harmonization increased tolerant Network Top1 from 4.7% to 20.3% and Network Top3 from 35.9% to 48.4%, but reduced tolerant broad-anatomy Top3 from 83.1% to 60.0%. Cohort-calibrated harmonization produced 20.3% Network Top1, 46.9% Network Top3 and 70.8% broad-anatomy Top3. Adapted preprocessing was not stable across independent EV-RNA and CSF-RNA representations. The locked baseline therefore remained the sole production route.

## Supplementary Figures

- **Figure S1.** Detailed internal validation across anatomical resolutions.
- **Figure S2.** Pairwise-rescue threshold and switch audit.
- **Figure S3.** AHBA cross-species validation and fine-resolution limitations.
- **Figure S4.** Glioma and liquid-biopsy transfer diagnostics.
- **Figure S5.** Exploratory adaptation sensitivity.

## Supplementary Tables

- **Table S1.** Datasets, sample sizes, label availability and analysis roles.
- **Table S2.** Internal validation metrics with confidence intervals.
- **Table S3.** External validation metrics.
- **Table S4.** Domain-adaptation sensitivity analysis.
- **Table S5.** Liquid-biopsy transfer stress tests.

## Submission checklist items

- Replace all public repository, archive, live-app, licence, author, affiliation, contact and funding placeholders.
- Create an immutable tagged release and archive it with a DOI.
- Ensure example data and tests run without private files.
- Upload the main figure as a separate high-resolution file and all supplementary figures/tables as one supplementary PDF or journal-supported package.
