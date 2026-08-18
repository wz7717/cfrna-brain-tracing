# Formal LOMO Exact F1 provenance

The current LOMO Exact endpoint is regenerated from the frozen prediction-level formal route.
Its macro denominator is the 104-label truth universe; Top1 labels outside that universe
remain false positives and do not create extra macro classes.

- Route: `top3_beam_local_top50_top100_zfusion_w0p25`
- Route family: `hybrid_projected_network_logcpm_exact`
- Source SHA-256: `EB5F10F01B122F68D09256EA6866DEAE2B439AABAD27E076181EC8760E7AAF36`
- Staged prediction SHA-256: `E90BB2ACF6C325530467A69252CBDF3B28C027EC41555FD03716527F5E3C101A`
- Top1 / Top3: `177/812`; `346/812`
- Prediction-only Top1 labels: `36c, Cla, MT`

## Derived summary

- macro_f1: `0.20341533397979938`
- sd_class_f1: `0.2159329084651328`
- median_class_f1: `0.15384615384615385`
- iqr_class_f1: `0.2857142857142857`
- weighted_f1: `0.21352686247014377`
- micro_f1: `0.21798029556650247`
- n_zero_f1_classes: `30`
- conditional_macro_f1_nonzero: `0.2858810099175559`

## Integer accounting

`sum(TP) = 177`; `sum(support) = 812`.
