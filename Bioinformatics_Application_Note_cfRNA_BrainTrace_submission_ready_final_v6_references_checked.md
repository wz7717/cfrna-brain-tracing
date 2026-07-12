cfRNA-BrainTrace: hierarchical brain-origin candidate ranking from RNA expression profiles with a primate transcriptomic atlas

Gene expression

Zhen Wang1, Li Zhang1, Zheng Wang2

1 Chinese Institute for Brain Research, Beijing, China
2 School of Psychological and Cognitive Sciences, Peking University, Beijing, China

To whom correspondence should be addressed: Li Zhang (zhangli@cibr.ac.cn) and Zheng Wang (zheng.wang@pku.edu.cn).

Abstract

Summary: cfRNA-BrainTrace is a Python and Streamlit application for hierarchical brain-origin candidate ranking from RNA expression profiles using a primate transcriptomic atlas. The locked production route uses projected variance-stabilized expression to generate a 10-class macaque Network Top3 beam, followed by logCPM-compatible reranking for resolution-group and exploratory exact-region candidates, while reporting marker coverage, entropy, score margins and scope warnings. In internal validation, Network Top3 accuracy reached 92.19% in leave-one-sample-out and 91.21% in leave-one-monkey-out evaluation across all 819 Network-evaluable samples, with lower accuracy at finer anatomical levels. External mapped-label and tumour/biofluid analyses define transfer limitations; cfRNA-BrainTrace is intended for reproducible coarse candidate ranking and resolution-limit auditing, not stand-alone clinical localization from unlabeled biofluid RNA.

Availability and implementation: cfRNA-BrainTrace is implemented in Python 3.11+ with command-line and Streamlit interfaces. Source code, README documentation, installation instructions, command-line and Streamlit entry points, unit tests, lightweight model artifacts, supplementary tables, validation/benchmark scripts and synthetic public example input/output files are available at https://github.com/wz7717/cfrna-brain-tracing. The manuscript-associated software release is v0.1.7 and is archived at Zenodo under the persistent concept DOI https://doi.org/10.5281/zenodo.20773674; the version DOI is assigned upon publication. A live Streamlit demonstration is available at https://brain-cfrna-tracing.streamlit.app/. The software is released under the MIT License.

Contact: zhangli@cibr.ac.cn; zheng.wang@pku.edu.cn

Supplementary information: Supplementary data are available at Bioinformatics online. The submitted supplementary file contains Supplementary Methods, Supplementary Results and embedded Supplementary Tables. The Markdown source, figure source data and synthetic public example input/output files are archived in the GitHub/Zenodo release.

1 Introduction

RNA expression profiles, including cell-free RNA profiles, can retain tissue- and cell-of-origin information (Vorperian et al., 2022), but within-brain tracing is limited by regional similarity, atlas granularity and domain shift. cfRNA-BrainTrace returns Network, resolution-group and exact-region candidates with confidence, coverage and scope diagnostics, emphasizing stable broad candidates when exact localization is not defensible.

2 System and methods

The locked production route was motivated by development analyses showing complementary behaviour of projected-VSD and logCPM-compatible expression spaces, and was then evaluated as a fixed three-tier route. Network denotes the 10-class macaque functional-anatomical source space used by the Bo2023 macaque transcriptomic reference (Bo et al., 2023): Cingulate gyrus; Frontal agranular motor areas; Hippocampal formation; Lateral Prefrontal Cortex; Occipital/Temporal; Operculum/Insula; Orbitomedial Prefrontal Cortex; Parietal and Parieto-occipital region; Subcortical; and Temporal. Projected-VSD denotes a projected variance-stabilized expression representation used for broad Network candidate generation, whereas downstream regional reranking is performed in a logCPM-compatible local expression space.

We used the Bo2023 macaque transcriptomic atlas because the task is hierarchical brain-origin ranking rather than a single human-atlas lookup. Bo2023 provides dense primate region-level RNA-seq, repeated animals and a curated Network-to-region hierarchy supporting leakage-free LOSO/LOMO validation, fold-local candidate generation and resolution-limit auditing. Human atlases, chiefly AHBA, were used only for mapped-label transfer tests.

Operationally, cfRNA-BrainTrace normalizes the uploaded expression table into logCPM- or logTPM-compatible space and aligns genes to the reference panel. A reference-fitted linear projector maps a query profile to Bo2023-like VSD space using stored gene-wise slope, intercept and clipping parameters from data/models/bo2023_reference_projector_linear_full.npz. This projection is single-sample compatible: it uses fixed reference-derived parameters and does not use target-cohort labels or cohort-level distribution information, reducing transductive-leakage risk. The projected representation is used only for Network Top3 beam generation. The top three Networks form the candidate beam, and regions outside this beam are excluded from downstream regional ranking. Within the retained beam, resolution groups and exact regions are reranked using logCPM-compatible local expression. The end-to-end workflow is illustrated in Fig. 1A.

For each query, expression values are represented as logCPM/logTPM-compatible values and, for the primary Network endpoint, projected into Bo2023-like VSD space. Let x denote the query vector, c_k the fold-local Network centroid and P(.) the fixed projection. Network candidates are ranked by s_k = corr(P(x), c_k). Pairwise rescue can only reorder Top1 within the retained Network Top3 beam: for high-confusion pairs (a,b), pair-specific training-fold markers define correlation scores q_a and q_b, and the larger score may move an already retained Network to rank 1. For downstream scoring, let B(x) be the retained beam and R(B) its candidate regions. Region candidates are reranked in logCPM-compatible space by z_r = lambda z{corr(x_G50, mu_r,G50)} + (1 - lambda) z{corr(x_G100, mu_r,G100)}, where mu_r is the fold-local centroid, G50/G100 are local discriminative panels, z{.} is within-candidate z-scoring and lambda=0.25 in the locked route. This value gives most weight to Top100 while retaining Top50 specificity and was fixed before formal P0 validation, not tuned on external datasets. Resolution-group scores use the same local evidence for groups not separable at high resolution.

Marker selection is performed only inside the appropriate training scope. The global Network panel ranks genes by a Fisher-like between-Network to within-Network variance ratio and retains the top 200 genes. Pairwise rescue uses 100 pair-specific discriminative genes for each high-confusion Network pair. Resolution-group reranking uses the top 200 fold-local or sample-local discriminative genes among candidate groups, whereas exploratory exact-region reranking fuses Top50 and Top100 local gene panels. Human/cross-species analyses use humanized signatures and report mapping coverage rather than treating unmapped macaque ENSMFAG-prefixed identifiers as human gene symbols.

3 Implementation

The Network route is implemented in core/network_tracing.py and the three-tier route in core/bo2023_region_tracing.py. The command-line and Streamlit interfaces call the same scoring core to avoid interface-specific divergence. Versioned artifacts in data/models/ store markers, centroids, anatomical dictionaries, projection parameters, route parameters and warning metadata. Unit tests cover Network scoring, region-resolution annotations, marker behaviour, upload metadata and VSD adaptation. Validation/export scripts and benchmark_runner.py support reruns of analyses, tables and figure artwork. Median post-reference-load inference time was 0.033 s per sample across 30 synthetic single-sample queries.

4 Validation

Validation was matched to available label resolution. The locked submission workflow uses projected-VSD representation for broad Network beam generation, followed by logCPM-compatible resolution-group and exploratory exact-region reranking within the retained Network beam. Exact-region output is therefore a candidate ranking rather than a deterministic localization endpoint.

The complete validation tested the three-tier architecture end to end: projected-VSD Network Top3 beam generation, logCPM-compatible resolution-group reranking and exact-region reranking. Network Top3 accuracy was 92.19% in leave-one-sample-out validation (755/819; 95% CI 90.14-93.83%) and 91.21% in leave-one-monkey-out validation (747/819; 95% CI 89.07-92.96%). Resolution-group Top3 accuracy was 72.36% among 814 LOSO reference-supported samples (95% CI 69.19-75.32%) and 69.09% among 812 LOMO reference-supported samples (95% CI 65.83-72.17%), whereas exact-region Top3 was 45.33% (95% CI 41.94-48.77%) and 42.36% (95% CI 39.01-45.79%) on the same denominators. Each internal endpoint exceeded its endpoint-specific uniform random baseline; for example, exact-region LOSO Top3 was 369/814 versus a 7.92% expectation (one-sided binomial p=4.43e-181). Exact McNemar/binomial discordant-pair tests showed no significant LOSO-vs-LOMO difference for Network Top3 (p=0.1686), a significant donor-level decrease for resolution-group Top3 (p=0.0389) and a borderline decrease for exact-region Top3 (p=0.0532). Five LOSO samples and seven LOMO samples lacked a truth-region reference after fold construction; they remained in Network evaluation but were not region-level errors or successes. The level-wise accuracy gradient supports Network Top3 as primary, resolution group as the main region-level output and exact region as exploratory local ranking. The validation summary is shown in Fig. 1B.

In Allen Human Brain Atlas (AHBA) mapped-label external validation (Hawrylycz et al., 2012), technical replicates were collapsed to independent tissue points and labels were harmonized using Network and region identifiers. Among network-qualified mapped-label points, the locked production route achieved Network Top1/Top3 accuracy of 73.99%/94.62% (`n=223`), resolution-group Top1/Top3 accuracy of 40.91%/62.50% (`n=88`) and exact-region Top1/Top3 accuracy of 23.86%/46.59% (`n=88`). These AHBA results should be interpreted as mapped-label transfer rather than direct anatomical equivalence, because human anatomical labels were harmonized to the macaque-derived hierarchy. In TCGA/BraTS glioma tissue RNA-seq with MRI-derived labels from the TCGA glioma MRI collections (Bakas et al., 2017), results support only coarse anatomical consistency because MRI truth labels are human atlas labels rather than Bo2023 macaque exact-region identifiers. After exclusion of one MRI-derived cerebellar out-of-scope case, Network Top3 was 40.63% (26/64) and broad-anatomy Top3 was 65.63% (42/64). The TCGA/BraTS results therefore support broad candidate consistency in tumour tissue but do not validate macaque Network-level or exact-region localization in human glioma. GSE189919 and other biofluid datasets without anatomical truth were treated as projection-feasibility or transfer stress tests rather than localization-accuracy validation.

5 Use and limitations

cfRNA-BrainTrace is intended for coarse brain-origin candidate ranking and resolution-aware assessment of whether a sample supports Network-level, resolution-group-level or only low-confidence output. It reports warnings for sparse profiles, low marker coverage, high entropy, low margins, out-of-scope anatomy and domain shift. Low-confidence samples should not be forced into exact-region predictions. cfRNA-BrainTrace is not intended for deterministic exact-region localization, clinical cfRNA localization without anatomical truth, or cerebellar/posterior-fossa localization outside the current reference space. Current evidence does not establish clinical liquid-biopsy localization.

The current reference and validation scope does not support cerebellar or posterior-fossa localization. The locked Bo2023-derived model was curated around available cerebral cortical, subcortical and related primate labels, and formal denominators include only labels represented in the training fold and mappable across datasets. Cerebellar and posterior-fossa structures were outside the validated label space, not hidden negative classes. Suspected posterior-fossa samples should be reported as out of scope. Future extension will require a cerebellum-inclusive reference, posterior-fossa label harmonization and patient-level anatomical truth.

The evidence should also be read against tissue-to-cfRNA domain shift. The model is trained and internally validated on brain tissue RNA-seq, whereas cfRNA profiles are sparse, mixture-like and affected by release, degradation, library preparation and non-brain background. Tissue LOSO/LOMO validation shows that the hierarchy is learnable in the reference domain, and AHBA/TCGA test transfer in human expression matrices. Biofluid datasets without patient-level anatomical truth are used only to audit gene coverage, projection feasibility and interpretable candidate distributions, preventing tissue-derived exact-region performance from being overgeneralized to clinical cfRNA localization.

Funding

This work was supported by the National Natural Science Foundation of China [grant number 32370682]; and the Capital Health Development Research Special Project [grant number Shou Fa 2026-1Q-1114].

Conflict of interest

None declared.

AI-assisted editing disclosure

The authors used large language model tools to assist with language editing, manuscript-format checking and code-documentation review. All scientific claims, analyses, code, validation results and final text were reviewed and approved by the authors, who take full responsibility for the content.

Data availability

The Bo2023 macaque transcriptomic atlas RNA-seq data are available through SRA accession PRJNA905082. The Allen Human Brain Atlas data are available from the Allen Brain Map portal. TCGA transcriptomic data are available from the NCI Genomic Data Commons, and the TCGA/BraTS MRI collections and segmentation resources are available from The Cancer Imaging Archive. GSE189919 is available from the NCBI Gene Expression Omnibus. Repository scripts, processed non-sensitive evaluation tables and figure source data are archived with the manuscript-associated v0.1.7 release under the project concept DOI https://doi.org/10.5281/zenodo.20773674; the release-specific DOI is assigned upon publication.

References

Bakas,S. et al. (2017) Advancing The Cancer Genome Atlas glioma MRI collections with expert segmentation labels and radiomic features. Sci. Data, 4, 170117.

Bo,T. et al. (2023) Brain-wide and cell-specific transcriptomic insights into MRI-derived cortical morphology in macaque monkeys. Nat. Commun., 14, 1499.

Hawrylycz,M.J. et al. (2012) An anatomically comprehensive atlas of the adult human brain transcriptome. Nature, 489, 391-399.

Vorperian,S.K. et al. (2022) Cell types of origin of the cell-free transcriptome. Nat. Biotechnol., 40, 855-861.

Figure 1. cfRNA-BrainTrace workflow and validation evidence. (A) cfRNA-BrainTrace takes an RNA expression matrix, accepts it through the command-line or Streamlit interface, projects the query into a macaque atlas-guided projected-VSD Network space to form a 10-class Network Top3 beam, and reranks resolution-group and exploratory exact-region candidates within the retained beam using logCPM-compatible local expression. Outputs are reported as hierarchical ranked candidates with confidence and diagnostics, including marker coverage, entropy, score margins and scope warnings. (B) Top3 validation performance for the locked route across internal leave-one-sample-out (LOSO) and leave-one-monkey-out (LOMO) validation, AHBA mapped-label validation and TCGA/BraTS coarse-consistency analysis. Internal and AHBA results show the expected anatomical-resolution gradient from Network to resolution group and exploratory exact region, whereas TCGA/BraTS supports only coarse anatomical consistency. Biofluid analyses lacking patient-level anatomical truth are not plotted as localization-accuracy results and are treated only as projection-feasibility or transfer-stress analyses.

Alt text: Two-panel scientific figure showing the cfRNA-BrainTrace workflow and validation summary. Panel A illustrates the analysis workflow: RNA expression matrix upload to cfRNA-BrainTrace, macaque brain atlas-guided projection, projected-VSD Network Top3 beam generation, logCPM-compatible local reranking and hierarchical candidate output with confidence and diagnostics. Panel B shows Top3 validation results. Internal LOSO Network, resolution-group and exact-region Top3 values are 92.19%, 72.36% and 45.33%; internal LOMO values are 91.21%, 69.09% and 42.36%; AHBA mapped-label values are 94.62%, 62.50% and 46.59%; TCGA/BraTS Network and broad-anatomy Top3 values are 40.63% and 65.63%. Biofluid datasets without patient-level anatomical truth are excluded from localization-accuracy bars.
