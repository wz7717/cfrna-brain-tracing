# Supplementary Material for cfRNA-BrainTrace

## Supplementary Methods

### S1 Reference atlas, gene-symbol audit and projector scope

The Bo2023 input resources were audited before manuscript use. Raw count, VSD and metadata matrices all contained 819 aligned samples. After symbol mapping, 21,668 common genes were available, and 199/200 locked Network-model genes were present in the common panel. Date-like Excel-mangled symbols were cleaned in the projector gene map and model artifact. The remaining suspicious identifiers were dominated by ENSMFAG fallback IDs or multi-ID aggregations and are reported as an interpretability limitation rather than a source-tracing failure.

### S2 Projector quality control

The projector was trained to map logCPM-compatible expression into Bo2023-like VSD space. Training-set sample-level reconstruction was strong (global Pearson 0.9961; median sample Pearson 0.9969; p10 sample Pearson 0.9940; MAE 0.2018; RMSE 0.2817). Gene-level agreement was weaker (median gene Pearson 0.6619; p10 gene Pearson 0.3665), so projected VSD was restricted to Network candidate generation.

### S3 Formal three-tier route

The production route first generates a Network Top3 beam in projected-VSD space. Candidate regions are then limited to those Networks. Resolution-group and exact-region candidates are reranked in logCPM-compatible local expression space using fold-local or route-local components. Accuracy is computed only when independent anatomical truth is available.

### S4 External validation policy

AHBA provides human normal brain RNA-seq with mapped labels and supports cross-species mapped-label validation. Exact-region metrics were calculated only for labels with stable Bo2023 mappings. TCGA/BraTS provides glioma tissue RNA-seq and corrected MRI-derived human labels; it supports only coarse anatomical consistency and is not a Bo2023 exact-region validation. GSE189919 was used only for projection feasibility because brain-origin accuracy labels were unavailable.

## Supplementary Results

### S1 Internal projected-VSD Network validation

In LOSO Network validation, projected VSD achieved Top1/Top3 accuracy of 0.5800/0.9158, compared with 0.5556/0.8828 for logCPM and 0.5311/0.8816 for native VSD. In LOMO validation, projected VSD achieved 0.5372/0.9133, compared with 0.4628/0.8510 for logCPM and 0.5043/0.8718 for native VSD.

### S2 Direct exact-region diagnostics

Direct exact-region scoring was weaker than Network scoring. In LOSO, projected-VSD exact Top1/Top3 was 0.1671/0.3563. In LOMO, projected-VSD exact Top1/Top3 was 0.1453/0.3116 and remained below native VSD Top3 of 0.3473. These diagnostics justify the hybrid downstream reranking route.

### S3 Formal three-tier validation

Formal LOSO hybrid validation yielded Network, group and exact Top3 accuracies of 0.9238, 0.7236 and 0.4533, respectively. Formal LOMO hybrid validation yielded 0.9121, 0.6909 and 0.4236. Exact-region performance remained clearly below Network performance, consistent with the intended interpretation hierarchy.

### S4 AHBA and special-label results

In AHBA mapped-label validation, the hybrid route achieved Network Top3 0.9442, Group Top3 0.6703 and Exact Top3 0.4286. Hybrid exact Top3 was higher than logCPM baseline (0.3077) and projected-VSD-only scoring (0.2967). For special labels, Insula, Caudate and Putamen all reached Network Top3 of 1.0000; Putamen exact Top3 was 0.8889, whereas Caudate exact Top1 was 0 and should not be described as an exact Top1 success.

### S5 TCGA/BraTS and GSE189919

In 65 TCGA/BraTS glioma tissue RNA-seq samples with MRI-derived labels, the hybrid route achieved Network Top3 0.4000 and Broad anatomy Top3 0.6462. These results support coarse anatomical consistency only. GSE189919 contained 51 samples and overlapped 15,622 of 21,668 projector genes (0.7210), supporting projection feasibility but not accuracy validation.

## Supplementary Tables

Table S1. Projector and gene-symbol audit.  
Table S2. Internal Network and formal three-tier validation.  
Table S3. Direct exact-region and local rerank diagnostics.  
Table S4. External AHBA, TCGA/BraTS and GSE189919 results.  
Table S5. Claim-boundary checklist for manuscript text.
