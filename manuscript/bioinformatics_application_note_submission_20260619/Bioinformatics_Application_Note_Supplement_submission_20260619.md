Supplementary Material for cfRNA-BrainTrace

Supplementary Methods

S1 Formal validation route

All validation analyses followed the same three-tier architecture used by the software. A query profile was first represented in logCPM- or logTPM-compatible expression space and then projected into Bo2023-like VSD space for Network scoring only. The top three Networks formed the candidate beam. Resolution-group and exact-region candidates were then reranked within that beam using local logCPM-compatible expression. Fold-local Network model construction differed between the LOSO and LOMO validation designs as documented in S7, but projected VSD remained restricted to Network beam generation in both. This design separates broad candidate generation from fine regional interpretation.

S1a Reference choice and cross-species interpretation

The choice of a macaque reference is methodological rather than a claim that macaque exact regions are directly interchangeable with human exact regions. Bo2023 provides a controlled primate RNA-seq reference with 819 samples, repeated donors and a dense set of cortical and subcortical labels that can be grouped into a stable Network, resolution-group and exact-region hierarchy. That hierarchy is essential for the tool: Network Top3 forms the primary beam, resolution groups provide the main region-level endpoint when local expression is not uniquely resolvable and exact-region output is deliberately exploratory. Human datasets are then used as transfer tests. AHBA supports mapped-label validation because its normal human brain labels can be harmonized to a subset of the macaque-derived hierarchy. Tumour and biofluid datasets support only coarser consistency, coverage or domain-transfer stress tests. This design avoids overfitting a sparse human validation label space while still testing whether the primate-derived signal transfers to real human expression matrices.

Important primate-to-human neuroanatomical differences remain. Human prefrontal association cortex occupies a different proportional and functional landscape from macaque prefrontal cortex, and language-dominant perisylvian systems do not have a one-to-one macaque homologue. These differences are why the model reports harmonized candidate rankings and resolution limits rather than asserting direct macaque-to-human exact-region equivalence.

S2 Internal validation design

Internal validation used the Bo2023 macaque brain reference. Two settings were evaluated. Leave-one-sample-out validation tested whether the route could recover the held-out sample label when other samples from the reference remained available. Leave-one-monkey-out validation held out all samples from one animal at a time and therefore tested donor-level generalization. Network metrics included all 819 samples because every held-out sample retained a supported Network label. Resolution-group and exact-region metrics required the truth region to remain represented in the corresponding training fold. Five LOSO samples and seven LOMO samples did not meet that region-reference requirement and were excluded only from region-level denominators. Top1, Top3 and median true-rank were used to distinguish exact calls from candidate-list recovery.

S2a Region-ranking score

The locked hybrid route is a constrained hierarchical ranker. The first level uses projected-VSD Network evidence only to define a small candidate beam. The second and third levels do not reuse projected VSD for exact localization; they instead rerank candidate regions with logCPM-compatible local expression. In fold f, the training samples define Network centroids c_k^(f), region centroids mu_r^(f) and a fold-local set of candidate regions R_f(B). Given query x_i, the Network beam is B_i = Top3_k corr(P_f(x_i), c_k^(f)). Pairwise rescue is restricted to high-confusion Network pairs already present in B_i. For a candidate pair (a,b), pair-specific gene set G_ab^(f) is selected in the training fold and the pairwise scores are q_a = corr(x_Gab, c_a,Gab^(f)) and q_b = corr(x_Gab, c_b,Gab^(f)). If q_a > q_b, Network a may be ranked ahead of b, and conversely for b; no Network outside the original Top3 beam can be introduced. The exact-region score for r in R_f(B_i) is a z-score fusion of correlations over two nested local gene panels, G50^(f) and G100^(f). The locked value lambda=0.25 weights the broader Top100 panel at 0.75 and the more selective Top50 panel at 0.25, reducing sensitivity to a very small marker set while retaining high-specificity local signal. This value was fixed in the hybrid route before the formal P0 validation refresh and was not optimized on AHBA, TCGA/BraTS or biofluid external datasets. This formulation separates the primary Network projection from downstream local expression evidence and prevents exact-region claims outside the candidate beam.

S3 External validation design

External validation was limited by the label resolution available in each dataset. AHBA human brain RNA-seq was used for mapped-label validation because its anatomical labels could be harmonized to Network, resolution-group and a subset of exact-region labels. TCGA/BraTS glioma tissue RNA-seq with MRI-derived labels was used only for coarse anatomical consistency because its truth labels are human imaging labels, not Bo2023 macaque exact-region identifiers. GSE189919 was used to test whether an external matrix could be projected into the model gene space; it was not used for accuracy estimation because patient-level anatomical truth was unavailable.

S4 Label harmonization and allowed conclusions

Label harmonization was performed only to the level supported by each dataset. AHBA anatomical labels were mapped to the macaque-derived Network, resolution-group and exact-region hierarchy where a supported mapping existed, and results are interpreted as mapped-label transfer rather than direct anatomical equivalence. TCGA/BraTS labels were derived from human MRI/tumour context and therefore support only coarse tumour-tissue anatomical consistency. Biofluid datasets without patient-level anatomical truth were not used for localization accuracy and are reported only as projection-feasibility or transfer stress tests.

For gene-symbol handling, macaque Ensembl fallback identifiers beginning with ENSMFAG were treated as valid macaque reference identifiers rather than spurious symbols. They were retained in internal macaque marker summaries but were not used as human gene symbols in human cfRNA or cross-species validation. Human-facing overlap analyses used the humanized signature file reports/bo2023_signature_genes_humanized_20260704.csv, in which 5324 of 8800 macaque signature rows (60.50%) mapped to a human gene symbol and 3476 remained unmapped ENSMFAG identifiers. Human validation and coverage statements therefore report the mapped humanized signature rather than counting unmapped macaque fallback IDs as human genes.

The 60.50% row-level humanization rate is a transfer constraint, not a failure of the macaque reference. It means that human-facing enrichment, coverage and cfRNA overlap analyses are conducted on the mapped subset and may underrepresent macaque-specific or poorly annotated genes. The formal AHBA validation is therefore interpreted at harmonized label levels and with explicit gene-coverage reporting, while unmapped ENSMFAG-prefixed rows are preserved only for macaque-reference transparency.

S5 Model artifacts and projected-VSD projector

The repository contains lightweight model artifacts in data/models/, including Network model files, region-resolution dictionaries, route metadata and the reference-fitted projector bo2023_reference_projector_linear_full.npz. The projector stores gene-wise slope, intercept and clipping parameters fitted from reference data. At inference time, a single query profile can be mapped into Bo2023-like projected-VSD space using these fixed parameters without target-cohort labels or target-cohort distribution fitting. Projected VSD is used only for Network Top3 beam generation; resolution-group and exploratory exact-region reranking use logCPM-compatible local expression within that beam.

The linear projection is an engineering compromise between computational efficiency and expression-space alignment. It is appropriate for a fixed, single-sample compatible transform, but it should not be interpreted as a complete model of cross-platform or cross-species expression distortion. Strong nonlinear batch effects, platform-specific gene-response curves or biological shifts not represented in the reference could reduce calibration; those cases should be handled through coverage, margin and entropy diagnostics and, where possible, dataset-specific validation rather than by assuming exact transfer.

S6 Example input/output

Synthetic public example files are provided in submission_ready_assets/example_io/ as part of the v0.1.6 public submission release. They document accepted input columns, Network ranked candidates, three-tier JSON output, resolution-group ranked candidates, exact-region ranked candidates and the local generation script. The example files are format and reproducibility aids only; they are not biological validation samples.

S6a Software availability

The manuscript-associated software release is v0.1.6. The public GitHub repository is available at https://github.com/wz7717/cfrna-brain-tracing, and the archived release is cited in the manuscript using the Zenodo version DOI https://doi.org/10.5281/zenodo.20780280; the project concept DOI is https://doi.org/10.5281/zenodo.20773674.

S7 Model development and locked evaluation timeline

The three-tier architecture was fixed before generating the final submission validation tables. External AHBA, TCGA/BraTS and biofluid analyses were not used to select the final route; they were used only for mapped-label transfer evaluation, coarse tumour-tissue consistency assessment and biofluid transfer stress testing. The LOSO implementation used locked Network genes and correlation-ranked projected-VSD Network Top3 candidates. The formal LOMO implementation rebuilt discriminative Network genes and a fold-local pairwise Top1 rescue within each training fold; this rescue could reorder Top1 but did not change the original Network Top3 candidate set. Both implementations used projected VSD only for the Network beam and logCPM-compatible local evidence for downstream reranking. The independent Network-only LOSO and LOMO analyses are reported separately as route-selection evidence.

S8 Figures and tables

Figure source data and synthetic example input/output files are archived in the GitHub/Zenodo release rather than treated as separate journal supplementary files. Supplementary Tables are embedded below so that the formal supplementary PDF can be submitted as a single file containing Supplementary Methods, Supplementary Results and the supporting tables. Together, these materials document both how the route was tested and what result each validation setting supports. Supplementary Figure S1 visualizes the three-tier hierarchy from 10 Networks to resolution groups and exploratory exact regions; the source file is manuscript/figures_publication/Supplementary_Figure_S1_three_tier_hierarchy.svg. Resolution-group names in Supplementary Figure S1 are derived from the curated Bo2023 hierarchy and the same bo2023_region_resolution_groups.json artifact used by the formal route.

S9 Evidence package provenance

The reviewer-response evidence package is organized into three tiers. P0 hard-evidence files contain the formal validation metrics, Wilson confidence intervals, random baselines, binomial tests, LOSO-vs-LOMO paired tests, confusion matrices, class-level F1, denominator audits, marker-methodology audits, Network anatomy tables, resolution-group hierarchy tables and AHBA mapping granularity summaries. P1 diagnostics contain dual-space consistency, error-cascade and confidence-diagnostic analyses. P2 completeness files contain development/comparator materials, including same-field tool positioning, simple ML baselines, test coverage, locked dependency files and random-seed registry. P2 comparator files are retained as development, comparator, sensitivity or diagnostic evidence; they are not merged into the formal three-tier hybrid validation endpoint.

Supplementary Results

S1 Internal route selection

The first internal analysis evaluated Network-level candidate generation. Projected VSD achieved Network Top1/Top3 of 58.00%/91.58% in LOSO and 53.72%/91.33% in LOMO, exceeding logCPM baseline and native VSD at Network Top3. Direct exact-region scoring was lower and less stable, especially in LOMO. These results motivated the use of projected VSD for broad Network beam generation and logCPM-compatible expression for downstream local reranking.

S2 Formal internal three-tier validation

In the complete LOSO validation, the formal route achieved Network Top1/Top3 of 58.24%/92.19% across all 819 samples. Resolution-group and exact-region Top3 were 72.36% and 45.33% among 814 reference-supported samples. In complete LOMO validation, Network Top3 was 91.21% across all 819 samples, while resolution-group and exact-region Top3 were 69.09% and 42.36% among 812 reference-supported samples. The 92.19% LOSO Network Top3 value uses all 819 Network-evaluable samples as the denominator. Region-level LOSO metrics use 814 reference-supported samples because five samples lacked a truth-region reference after fold construction. Region-level LOMO metrics use 812 reference-supported samples because seven samples lacked a truth-region reference after fold construction. The earlier LOSO Network value of 92.38% was conditional on the 814 region-evaluable samples and is retained only as a legacy denominator inconsistency, not as the submission result. Median true-rank increased from Network to exact-region levels, consistent with decreasing anatomical certainty at finer resolution. Resolution group is therefore the preferred region-level endpoint, while exact-region output is retained as a candidate ranking.

Wilson confidence intervals and statistical tests were generated for the formal P0 evidence package. Network Top3 was 92.19% in LOSO (755/819; 95% CI 90.14-93.83%) and 91.21% in LOMO (747/819; 95% CI 89.07-92.96%). Resolution-group Top3 was 72.36% in LOSO (589/814; 95% CI 69.19-75.32%) and 69.09% in LOMO (561/812; 95% CI 65.83-72.17%). Exact-region Top3 was 45.33% in LOSO (369/814; 95% CI 41.94-48.77%) and 42.36% in LOMO (344/812; 95% CI 39.01-45.79%). Each endpoint exceeded its uniform random baseline by one-sided binomial testing, including exact-region LOSO Top3 versus 7.92% uniform random expectation (p=4.43e-181). Paired LOSO-vs-LOMO tests showed no significant Network Top3 decrease (p=0.1686), a significant resolution-group decrease (p=0.0389) and a borderline exact-region decrease (p=0.0532).

S3 External validation results

In AHBA, the formal three-tier hybrid route achieved Network Top1/Top3 of 74.68% (95% CI 68.73-79.83%)/94.42% (90.69-96.71%), resolution-group Top1/Top3 of 36.26% (27.13-46.51%)/67.03% (56.86-75.83%) and exact-region Top1/Top3 of 24.18% (16.54-33.90%)/42.86% (33.18-53.11%). Hybrid exact Top3 exceeded logCPM baseline (30.77%) and projected-VSD-only scoring (29.67%). These AHBA values are the main external validation result because the current manuscript uses AHBA as mapped-label human brain transfer, not as a historical route-development trace. In TCGA/BraTS, hybrid Network Top3 was 40.00% and broad-anatomy Top3 was 64.62%, supporting only coarse anatomical consistency. Against a nominal 30% Network Top3 reference level, the TCGA/BraTS Network Top3 value of 26/65 has a one-sided normal-approximation p value of approximately 0.039; the more conservative exact binomial test gives p=0.0548. Against the endpoint-specific weighted random baseline of 20.28%, the exact one-sided binomial p value is 2.17e-4. These values support only coarse tumour-tissue consistency and do not convert TCGA/BraTS into an exact localization validation set. GSE189919 overlapped 15,622/21,668 projector genes, corresponding to 72.10% gene-space coverage, supporting projection feasibility rather than source-localization accuracy.

S4 Dual-space and diagnostic analyses

Projected VSD and logCPM are not independent competing endpoints. Projected VSD defines the broad Network Top3 candidate beam, and logCPM then performs local reranking inside that beam. In the P1 diagnostic package, LOSO exact-region Top1 candidates were inside the retained Network beam in 100.00% of cases, LOSO exact-region Top3 candidates were all inside the beam in 100.00% of cases and no Top3 candidate occurred outside the retained beam. Conditional on a correct Network beam, LOSO exact-region Top3 was 49.07% and resolution-group Top3 was 78.32%; after a missed Network beam, LOSO exact-region recovery was 0.00%. These results support the intended cascade interpretation: Network misses are near-hard failures for downstream recovery, whereas remaining downstream misses after a correct beam reflect local anatomical ambiguity.

Diagnostic outputs were evaluated as confidence signals rather than as independent validation endpoints. Full Network scores provided top1_margin, top3_beam_margin, normalized score entropy and marker-overlap counts. For LOSO exact-region Top3, the Top3 beam margin had AUC 0.613 for downstream hit/miss discrimination, whereas marker overlap was not testable in internal Bo2023 or AHBA validation because the relevant marker panels were fully covered. An AUC of 0.613 is only modestly above the random-discrimination value of 0.5, so score margin should not be used as a binary accept/reject decision rule; it is a weak confidence aid to be interpreted together with anatomical level, marker coverage, entropy and dataset scope. External biofluid matrices without anatomical truth therefore support only coverage and projection-feasibility checks. The GSE189919 projector overlap result is used in this limited sense: it indirectly tests whether a sparse external profile can enter the gene space, but it does not test cfRNA degradation robustness or localization accuracy.

Reported p values are per-endpoint and were not adjusted for multiple comparisons because the endpoints address distinct anatomical-resolution hypotheses: Network candidate recovery, resolution-group recovery, exploratory exact-region recovery and coarse external consistency.

S5 Comparator, sensitivity and reproducibility materials

Simple ML baselines were run as leave-one-monkey-out Network classifiers on Bo2023 VSD expression with training-fold feature selection. Top3 accuracies were 80.46% for 5-nearest-neighbour cosine classification, 83.27% for nearest centroid and 90.96% for random forest. These baselines are retained as comparator evidence only: they do not return resolution groups, exact-region candidate rankings, resolution-limit diagnostics or the formal three-tier route used by the software. Same-field comparison found no strict like-for-like public tool. The closest named comparators are CIBERSORTx and related cell-type deconvolution tools, TissueEnrich and tissue-expression enrichment tools, broad cell-free transcriptome tissue-of-origin studies, and Allen/AHBA atlas-query workflows. cfRNA-BrainTrace should therefore be positioned as a resolution-aware hierarchical candidate-ranking and audit tool rather than a generic tissue deconvolution method.

The baseline comparison supports the formal route rather than replacing it. Random forest approaches the formal LOMO Network Top3 result at the broad Network level, but it does not provide the required hierarchical outputs, fold-local resolution-group reranking, exploratory exact-region candidate lists, score-margin/entropy diagnostics or explicit abstention boundaries. The formal hybrid route is therefore preferred for the manuscript because it is a validated end-to-end anatomical ranking workflow, whereas ML baselines are restricted comparator checks for broad Network classification.

Formal LOSO sensitivity analyses reran the same P0 three-tier pipeline with alternative local exact-region fusion weights. The locked route used lambda=0.25, which gave the highest exact-region Top1 accuracy among the evaluated values; lambda=0.50 gave a numerically higher Top3 value, while lambda=0.75 was slightly lower than lambda=0.50. Because the manuscript primary route had been locked before the formal validation refresh and because all three values gave similar median true-rank, lambda=0.25 was retained as the prespecified formal route rather than reselected by post hoc Top3 optimization.

The engineering reproducibility package includes a pinned requirements-lock.txt, an environment.yml, a random-seed registry and coverage artifacts. These files address software reproducibility and reviewer traceability but do not alter the formal validation route.

S6 Cerebellum/posterior-fossa scope and cfRNA domain transfer

We intentionally avoid inferring cerebellar or posterior-fossa origin from the present model. The reason is not that these locations are biologically irrelevant, but that the current locked evidence chain lacks a validated cerebellum-inclusive label space. A classifier cannot be judged on a label class that is absent from the training reference, absent from the candidate hierarchy or not represented in the external truth mapping. Including such labels only as textual possibilities would create an apparent localization capability without an evaluable denominator. The correct treatment is therefore abstention/out-of-scope reporting. Future work should add cerebellar transcriptomic reference profiles, posterior-fossa label harmonization, MRI-linked posterior-fossa truth and a prospective validation design before any cerebellar or posterior-fossa performance claim is made.

Several layers of transfer separate the training reference from the intended cfRNA use case. First, the reference is tissue RNA-seq, while cfRNA reflects extracellular RNA abundance after cell release, degradation and clearance. Second, cfRNA contains contributions from multiple tissues and blood components, so the brain signal may be diluted. Third, external cfRNA studies often lack patient-level anatomical truth, making apparent localization impossible to validate directly. We therefore separate evidence types: internal Bo2023 LOSO/LOMO measures reference-domain traceability; AHBA measures mapped-label transfer to normal human brain RNA-seq; TCGA/BraTS measures coarse consistency in tumour tissue with MRI-derived labels; and biofluid datasets measure gene coverage and projection feasibility only. The required future validation is a prospective or curated cfRNA cohort with matched imaging or surgical anatomical truth, reported at the same Network/resolution-group/exact-region levels and with abstention for unsupported labels.

Supplementary Tables

Table S1. Internal validation design

| Validation | Data | Held-out unit | Route tested | Reported endpoints | Denominator policy |
|---|---|---|---|---|---|
| Network LOSO | Bo2023 macaque brain RNA-seq | Single sample | Projected VSD vs logCPM/native VSD Network scoring | Network Top1, Network Top3, median true-rank | All 819 samples |
| Network LOMO | Bo2023 macaque brain RNA-seq | One monkey | Projected VSD vs logCPM/native VSD Network scoring | Network Top1, Network Top3, median true-rank | All 819 samples |
| Formal three-tier LOSO | Bo2023 macaque brain RNA-seq | Single sample | Projected-VSD Network beam plus logCPM group/exact rerank | Network, resolution-group and exact-region Top1/Top3 | Network n=819; group/exact n=814 |
| Formal three-tier LOMO | Bo2023 macaque brain RNA-seq | One monkey | Fold-local Network beam plus logCPM group/exact rerank | Network, resolution-group and exact-region Top1/Top3 | Network n=819; group/exact n=812 |

Table S1a. Ten Network labels and anatomy mapping

| Network | Regions | Resolution groups | Training samples | Classical neuroanatomy mapping |
|---|---:|---:|---:|---|
| Cingulate gyrus | 10 | 8 | 39 | Medial cingulate cortex; limbic/medial association cortex. |
| Frontal (agranular frontal motor areas) | 7 | 3 | 65 | Agranular frontal and premotor cortex. |
| Hippocampal formation | 1 | 1 | 8 | Hippocampal/parahippocampal formation. |
| Lateral Prefrontal Cortex | 14 | 4 | 108 | Dorsolateral and ventrolateral prefrontal association cortex. |
| Occipital/Temporal | 9 | 5 | 69 | Visual occipital cortex and occipito-temporal visual association cortex. |
| Operculum/Insula | 11 | 2 | 81 | Opercular and insular cortex. |
| Orbitomedial Prefrontal Cortex (OMPFC) | 13 | 4 | 103 | Orbitofrontal and medial prefrontal cortex. |
| Parietal, and Parieto-occipital region | 11 | 4 | 99 | Posterior parietal and parieto-occipital association cortex. |
| Subcortical | 9 | 7 | 54 | Basal ganglia, amygdala, thalamic and related subcortical labels in this reference. |
| Temporal | 27 | 12 | 193 | Temporal association and auditory-related cortex. |

Table S2. Internal validation results

| Dataset | Route | Endpoint | Evaluated n | Top1 | Top3 | Median true-rank | Interpretation |
|---|---|---|---:|---|---|---|---|
| Bo2023 LOSO | Projected VSD Network | Network | 819 | 58.00% | 91.58% | 1.0 | Supports projected-VSD Network beam |
| Bo2023 LOMO | Projected VSD Network | Network | 819 | 53.72% | 91.33% | 1.0 | Supports donor-level Network beam |
| Bo2023 LOSO | Formal hybrid | Network | 819 | 58.24% | 92.19% | 1.0 | Primary endpoint across all Network-evaluable samples |
| Bo2023 LOSO | Formal hybrid | Resolution group | 814 | 44.47% | 72.36% | 2.0 | Main region-level endpoint among reference-supported samples |
| Bo2023 LOSO | Formal hybrid | Exact region | 814 | 22.48% | 45.33% | 4.0 | Exploratory ranking among reference-supported samples |
| Bo2023 LOMO | Formal hybrid | Network | 819 | 57.75% | 91.21% | 1.0 | Cross-monkey support |
| Bo2023 LOMO | Formal hybrid | Resolution group | 812 | 41.38% | 69.09% | 2.0 | Main region-level endpoint among reference-supported samples |
| Bo2023 LOMO | Formal hybrid | Exact region | 812 | 22.17% | 42.36% | 5.0 | Exploratory ranking among reference-supported samples |

Table S3. External validation design

| Dataset | Sample type | n | Truth label type | Allowed conclusion |
|---|---|---|---|---|
| AHBA | Human normal brain RNA-seq | 242 total; 233 supported; 91 exact-evaluable | Mapped anatomical labels | Cross-species mapped-label validation |
| TCGA/BraTS | Glioma tissue RNA-seq with MRI-derived labels | 65 patients | Human imaging labels | Coarse anatomical consistency only |
| GSE189919 | External count matrix | 51 samples | No patient-level anatomical truth | Projection feasibility only |

Table S4. External validation results

| Dataset | Route | Endpoint | Top1 | Top3 | Conclusion |
|---|---|---|---|---|---|
| AHBA | Hybrid | Network | 74.68% | 94.42% | Strong mapped-label Network support |
| AHBA | Hybrid | Resolution group | 36.26% | 67.03% | Moderate group-level support |
| AHBA | Hybrid | Exact region | 24.18% | 42.86% | Only exact-evaluable mapped labels |
| AHBA | logCPM baseline | Exact region | 17.58% | 30.77% | Below hybrid exact Top3 |
| AHBA | Projected VSD only | Exact region | 10.99% | 29.67% | Below hybrid exact Top3 |
| TCGA/BraTS | Hybrid | Network | 15.38% | 40.00% | Coarse consistency only |
| TCGA/BraTS | Hybrid | Broad anatomy | 13.85% | 64.62% | Coarse consistency only |
| GSE189919 | Projection feasibility | Projector gene overlap | 15622/21668 | 72.10% | No accuracy claim |

Table S4a. AHBA formal hybrid 95% confidence intervals and random baselines

| Endpoint | Top1, % (95% CI) | Top3, % (95% CI) | Uniform Top3 baseline | Weighted Top3 baseline |
|---|---:|---:|---:|---:|
| Network, n=233 | 74.68 (68.73-79.83) | 94.42 (90.69-96.71) | 37.5% | 56.9% |
| Resolution group, n=91 | 36.26 (27.13-46.51) | 67.03 (56.86-75.83) | 20.0% | 30.2% |
| Exact region, n=91 | 24.18 (16.54-33.90) | 42.86 (33.18-53.11) | 30.0% | 32.1% |

Table S4b. Internal weighted random baselines

| Endpoint | Uniform Top3 baseline | Weighted Top3 baseline | Observed Top3 |
|---|---:|---:|---:|
| Formal Network LOSO | 30.0% | 38.4% | 92.19% |
| Formal Network LOMO | 30.0% | 38.4% | 91.21% |
| Formal resolution-group LOSO | 22.0% | 6.3% | 72.36% |
| Formal resolution-group LOMO | 20.4% | 4.4% | 69.09% |
| Formal exact-region LOSO | 7.9% | 3.2% | 45.33% |
| Formal exact-region LOMO | 7.9% | 3.2% | 42.36% |

Table S5. Figure and table index

| Item | Content | Use in manuscript |
|---|---|---|
| Figure 1 | Two-panel workflow and validation summary; Panel B numeric source in Figure1_validation_summary.csv | Main Application Note figure |
| Table S1 | Internal validation design | Documents LOSO/LOMO operations |
| Table S1a | Ten Network labels and anatomy mapping | Lists all 10 Network classes, region counts, resolution-group counts and classical neuroanatomy mapping |
| Table S2 | Internal validation results | Reports Network, group and exact-region metrics |
| Table S3 | External validation design | Documents label support and allowed conclusions |
| Table S4 | External validation results | Reports AHBA, TCGA/BraTS and GSE189919 outcomes |
| Table S4a | AHBA formal hybrid 95% confidence intervals and random baselines | Documents uncertainty and baseline context for the main AHBA external validation |
| Table S4b | Internal weighted random baselines | Documents endpoint-specific uniform and weighted random baselines for formal internal validation |
| Table S5 | Figure/table index | Indexes manuscript figure and supplementary tables |
| Table S6 | Claim boundaries | Defines unsupported interpretations and permitted claims |
| Table S7 | Marker selection methodology | Documents marker selection stages, feature counts and fold policy |
| Table S8 | Network-level class F1 | Reports class-level precision, recall and F1 for Network endpoints |
| Table S9 | Largest Network Top1 confusion pairs | Summarizes major Network-level error modes and points to repository supplement matrices |
| Table S10 | Named same-field comparator families | Positions cfRNA-BrainTrace relative to deconvolution, enrichment and atlas-query comparators |
| Table S11 | Local exact-region fusion-weight sensitivity | Reports available archived lambda sensitivity values for exact-region local fusion |

Table S6. Claim boundaries

| Avoid | Use |
|---|---|
| Projected VSD creates a new Bo2023 atlas. | Only query profiles are projected for Network beam generation. |
| Projected VSD is best for exact-region inference. | Projected VSD supports Network beam; logCPM supports downstream reranking. |
| TCGA/BraTS validates Bo2023 exact regions. | TCGA/BraTS supports coarse anatomical consistency. |
| GSE189919 validates accuracy. | GSE189919 verifies projection feasibility only. |
| Exact Top1 is deterministic localization. | Exact outputs are exploratory local candidate rankings. |

Table S7. Marker selection methodology

| Component | Stage | Selection rule | Features | Fold policy |
|---|---|---|---:|---|
| Network global marker panel | Network Top3 beam generation | Fisher-like between-Network / within-Network variance-ratio ranking on Bo2023 training data | 200 | Rebuilt fold-locally in validation; locked full-reference panel in production |
| Network pairwise rescue markers | Top1 reordering within retained Network Top3 beam | High-confusion Network-pair discriminative genes | 100 per pair | Rebuilt fold-locally in LOMO formal validation; constrained to original beam |
| Resolution-group local genes | Group reranking inside Network Top3 beam | Local discriminative ranking among candidate groups | 200 | Fold-local/sample-local; held-out truth excluded |
| Exact-region local genes | Exploratory exact-region reranking | Nested local Top50/Top100 correlation z-score fusion | 50+100 fused | Fold-local/sample-local; exact endpoint exploratory |
| Development marker annotation route | Diagnostic support only | Stable markers with min_consistency=0.75, min_effect=0.5 and min_markers_for_support=8 | Top30 per region | Not adopted as main reranking route |

Table S8. Network-level class F1

| Network | LOSO support | LOSO precision | LOSO recall | LOSO F1 | LOMO precision | LOMO recall | LOMO F1 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Cingulate gyrus | 39 | 0.51 | 0.59 | 0.55 | 0.56 | 0.62 | 0.59 |
| Frontal (agranular frontal motor areas) | 65 | 0.59 | 0.54 | 0.56 | 0.65 | 0.68 | 0.66 |
| Hippocampal formation | 8 | 0.57 | 1.00 | 0.73 | 0.67 | 1.00 | 0.80 |
| Lateral Prefrontal Cortex | 108 | 0.57 | 0.47 | 0.52 | 0.67 | 0.43 | 0.52 |
| Occipital/Temporal | 69 | 0.71 | 0.71 | 0.71 | 0.73 | 0.70 | 0.71 |
| Operculum/Insula | 81 | 0.42 | 0.57 | 0.48 | 0.39 | 0.65 | 0.49 |
| Orbitomedial Prefrontal Cortex (OMPFC) | 103 | 0.66 | 0.59 | 0.62 | 0.61 | 0.50 | 0.55 |
| Parietal, and Parieto-occipital region | 99 | 0.43 | 0.78 | 0.56 | 0.43 | 0.74 | 0.54 |
| Subcortical | 54 | 1.00 | 0.78 | 0.88 | 1.00 | 0.78 | 0.88 |
| Temporal | 193 | 0.71 | 0.44 | 0.54 | 0.65 | 0.44 | 0.52 |

Table S9. Largest Network Top1 confusion pairs

| Validation | True Network | Predicted Network | Count |
|---|---|---|---:|
| LOSO | Temporal | Operculum/Insula | 47 |
| LOSO | Temporal | Parietal, and Parieto-occipital region | 46 |
| LOSO | Orbitomedial Prefrontal Cortex (OMPFC) | Lateral Prefrontal Cortex | 21 |
| LOSO | Lateral Prefrontal Cortex | Frontal (agranular frontal motor areas) | 21 |
| LOSO | Operculum/Insula | Temporal | 15 |
| LOSO | Lateral Prefrontal Cortex | Parietal, and Parieto-occipital region | 15 |
| LOSO | Occipital/Temporal | Parietal, and Parieto-occipital region | 14 |
| LOSO | Temporal | Occipital/Temporal | 13 |
| LOMO | Temporal | Parietal, and Parieto-occipital region | 49 |
| LOMO | Temporal | Operculum/Insula | 49 |
| LOMO | Orbitomedial Prefrontal Cortex (OMPFC) | Operculum/Insula | 20 |
| LOMO | Lateral Prefrontal Cortex | Frontal (agranular frontal motor areas) | 20 |
| LOMO | Lateral Prefrontal Cortex | Parietal, and Parieto-occipital region | 17 |
| LOMO | Lateral Prefrontal Cortex | Orbitomedial Prefrontal Cortex (OMPFC) | 15 |
| LOMO | Orbitomedial Prefrontal Cortex (OMPFC) | Lateral Prefrontal Cortex | 14 |
| LOMO | Operculum/Insula | Temporal | 14 |

Full machine-readable confusion matrices and class-level F1 tables for Network, resolution-group and exact-region endpoints are available in the repository supplement.

Table S10. Named same-field comparator families

| Tool or family | Primary goal | Difference from cfRNA-BrainTrace |
|---|---|---|
| CIBERSORTx and related cell-type deconvolution tools | Estimate cell-type fractions from bulk expression | Cell-type composition, not brain-region hierarchical source ranking |
| TissueEnrich and tissue-expression enrichment tools | Test tissue enrichment from gene sets or expression signatures | Broad tissue-level inference, not Network to resolution-group to exact-region candidate ranking |
| Cell-free transcriptome tissue-of-origin studies | Infer broad tissue or cell-type contribution to cfRNA | Important cfRNA context but usually not within-brain atlas localization |
| Allen/AHBA atlas-query workflows | Map genes or samples to human brain atlas annotations | Atlas-query context; not packaged as a validated app/CLI hierarchical ranker |
| cfRNA-BrainTrace | Hierarchical brain-origin candidate ranking and resolution-limit auditing | Explicit Top3 beam, three-tier output, confidence diagnostics and shared Streamlit/CLI scoring core |

Table S11. Local exact-region fusion-weight sensitivity

| Route variant | Evaluated n | Exact-region Top1 | Exact-region Top3 | Median true-rank | Interpretation |
|---|---:|---:|---:|---:|---|
| lambda=0.25 | 814 | 22.48% (183/814) | 45.33% (369/814) | 4.0 | Locked formal route; highest Top1 among evaluated formal LOSO weights |
| lambda=0.50 | 814 | 22.36% (182/814) | 46.07% (375/814) | 4.0 | Slightly higher Top3 but lower Top1 than the locked route |
| lambda=0.75 | 814 | 22.24% (181/814) | 45.45% (370/814) | 4.0 | Similar median rank; lower Top1 than the locked route |
