| Component | Implementation | Principal output |
| --- | --- | --- |
| Input | Gene-symbol expression table; optional metadata | Marker overlap and non-zero coverage |
| Scoring | 200-gene Pearson correlation plus pairwise rescue | Ranked Network Top1/Top3 and rescue status |
| Hierarchy | Lobe, broad anatomy, Network and exact region | Resolution-specific candidates and out-of-scope flags |
| Diagnostics | Correlation, display probability, margin and entropy | Confidence and transfer-quality audit |
| Interfaces | Streamlit web app and command-line interface | Interactive review and reproducible batch export |
| Validation policy | Accuracy only with independent anatomical truth | Unlabelled biofluids reported as stress tests |
