# 本线程测试结果与论文初稿交接指导

生成日期：2026-06-19

用途：交给另一个线程撰写论文初稿。本文只总结本线程中已经落盘、有文件证据支持的测试结果；没有结果文件支持的内容不写成结论。

## 0. 总结结论

当前可写入论文初稿的核心路线是：

```text
projected VSD Network Top3 beam
  -> logCPM resolution group rerank
  -> logCPM local exact rerank
```

精确解释：

- 不是把整个 Bo2023 猕猴脑图谱重新投影成一个新 atlas。
- Network 层训练参考仍使用 Bo2023 原版 VSD/batch-removed 参考矩阵。
- LOSO/LOMO 中被留出的样本或外部 query 会从 logCPM/logTPM 通过 projector 映射到 Bo2023-like VSD 空间，用于 Network Top3 candidate beam。
- resolution group 与 exact region 层回到 logCPM 局部表达空间进行重排序。

论文中应把 `Network Top3` 作为主 endpoint，把 `resolution group` 作为 region-level 主要报告粒度，把 `exact region` 作为局部候选排序。不要把 exact Top1 写成确定性结论。

## 1. 数据与 projector 质量审计

证据文件：

- `results/bo2023_reference_projection_20260616_cleaned_symbols/data_audit_summary.json`
- `results/bo2023_reference_projection_20260616_cleaned_symbols/projector_qc_summary.json`
- `results/bo2023_reference_projection_20260616_cleaned_symbols/gene_symbol_audit_summary.json`

### 1.1 Bo2023 counts/VSD/metadata 对齐

已确认：

- raw count genes: `28415`
- raw VSD genes: `23605`
- count gene symbols: `26478`
- VSD gene symbols: `21668`
- count samples: `819`
- VSD samples: `819`
- metadata samples: `819`
- common samples: `819`
- common genes after symbol mapping: `21668`
- locked Network model genes: `200`
- locked Network model genes in common panel: `199`

重要边界：

- `missing_locked_model_genes` 中曾出现 `2021-09-06 00:00:00`，这是 Excel 日期样式 gene symbol 问题的一部分。
- 后续 cleaned gene map 与模型 artifact 已修复日期样式 gene symbol；不要在论文中把日期样式 symbol 当作真实基因。

### 1.2 Projector 拟合质量

projector 训练集拟合指标：

| 指标 | 数值 |
|---|---:|
| n_genes | 21668 |
| n_samples | 819 |
| global Pearson | 0.9961 |
| MAE | 0.2018 |
| RMSE | 0.2817 |
| median sample Pearson | 0.9969 |
| p10 sample Pearson | 0.9940 |
| median gene Pearson | 0.6619 |
| p10 gene Pearson | 0.3665 |
| median gene R2 | 0.4402 |
| p10 gene R2 | 0.1343 |
| fallback genes | 0 |

可写结论：

- projector 能很好地重建样本层面的 VSD 表达模式，适合用于 Network-level candidate generation。
- gene-level 相关性分布比 sample-level 弱，因此不应声称 projector 在每个基因上都精确重建 VSD。

不要写：

- “projector 完美恢复每个基因的 VSD。”
- “projector 生成了新的 Bo2023 atlas。”

### 1.3 Gene symbol 审计

gene symbol 审计结果：

- projector genes: `21668`
- ok symbols: `16169`
- suspicious symbols: `5499`
- date-like symbols: `0`
- Ensembl fallback symbols: `5499`
- suspicious in locked Network model: `19`
- date-like in locked Network model: `0`

解释：

- 日期样式 Excel-mangled symbols 已清理。
- 剩余 suspicious symbols 主要是 `ENSMFAG...` fallback IDs 或 multi-ID aggregation，不是日期转换 bug。
- fallback IDs 可能降低外部生物学可解释性，但不等于样本无法溯源；论文中应报告或说明 fallback-ID fraction。

## 2. 外部投影可行性测试

证据文件：

- `results/bo2023_reference_projection_20260616_cleaned_symbols/external_projected_vsd_GSE189919_summary.json`

GSE189919 projection 可行性：

- input genes: `59453`
- projector genes: `21668`
- overlap projector genes: `15622`
- overlap fraction: `0.7210`
- samples: `51`

可写结论：

- 外部 raw count matrix 可以映射到 projector gene space，说明 reference projection 在基因覆盖层面可执行。

不要写：

- GSE189919 在本线程中没有可用于 brain-origin accuracy 的标签结果；不要把它写成准确率验证。

## 3. 内部 Network 层测试

证据文件：

- `results/bo2023_reference_projection_20260616_cleaned_symbols/bo2023_projected_vsd_loso_summary.json`
- `results/bo2023_reference_projection_20260616_cleaned_symbols/bo2023_projected_vsd_lomo_summary.json`

### 3.1 LOSO Network

测试设计：fold-local SaleemNetworks LOSO，使用 locked Network-model genes，`n=819`。

| Route | Network Top1 | Network Top3 | Median true rank |
|---|---:|---:|---:|
| logCPM baseline | 0.5556 | 0.8828 | 1.0 |
| native VSD | 0.5311 | 0.8816 | 1.0 |
| projected VSD | 0.5800 | 0.9158 | 1.0 |

正向结果：

- projected VSD 在 LOSO Network Top1/Top3 均高于 logCPM baseline 和 native VSD。
- 该结果支持把 projected VSD 用作 Network Top3 beam。

### 3.2 LOMO Network

测试设计：strict leave-one-monkey-out SaleemNetworks validation，`n=819`。

| Route | Network Top1 | Network Top3 | Median true rank |
|---|---:|---:|---:|
| logCPM baseline | 0.4628 | 0.8510 | 2.0 |
| native VSD | 0.5043 | 0.8718 | 1.0 |
| projected VSD | 0.5372 | 0.9133 | 1.0 |

正向结果：

- projected VSD 在跨猴泛化中 Network Top3 达到 `0.9133`，明显高于 logCPM baseline 和 native VSD。
- 该结果是选择 projected VSD 作为 Network candidate beam 的最直接内部证据。

## 4. Direct exact-region 诊断测试

证据文件：

- `results/bo2023_reference_projection_20260616_cleaned_symbols/bo2023_projected_vsd_exact_region_loso_summary.json`
- `results/bo2023_reference_projection_20260616_cleaned_symbols/bo2023_projected_vsd_exact_region_lomo_summary.json`

### 4.1 LOSO direct exact

| Route | Valid n | Exact Top1 | Exact Top3 | Median true rank |
|---|---:|---:|---:|---:|
| logCPM baseline | 814 | 0.1376 | 0.2604 | 13.0 |
| native VSD | 814 | 0.1671 | 0.3403 | 7.0 |
| projected VSD | 814 | 0.1671 | 0.3563 | 6.0 |

正向结果：

- direct projected VSD exact Top3 在 LOSO 中高于 logCPM baseline 和 native VSD。

限制：

- 即使最好也只有 exact Top3 `0.3563`，不足以作为主线 exact-region 终点。

### 4.2 LOMO direct exact

| Route | Valid n | Exact Top1 | Exact Top3 | Median true rank |
|---|---:|---:|---:|---:|
| logCPM baseline | 812 | 0.1392 | 0.2537 | 13.0 |
| native VSD | 812 | 0.1773 | 0.3473 | 7.0 |
| projected VSD | 812 | 0.1453 | 0.3116 | 8.0 |

关键结论：

- direct projected VSD 在 LOMO exact-region 层低于 native VSD。
- 因此不能写 “projected VSD 在 exact-region 层优于所有路线”。
- direct exact scoring 应作为诊断基线，不应作为最终主线。

## 5. Region-local rerank 测试

证据文件：

- `results/bo2023_reference_projection_20260616_cleaned_symbols/bo2023_projected_vsd_region_local_rerank_summary.json`
- `results/bo2023_reference_projection_20260616_cleaned_symbols/bo2023_projected_vsd_region_local_rerank_loso_lomo_summary.json`
- `results/bo2023_reference_projection_20260616_cleaned_symbols/region_local_rerank_loso_hybrid/bo2023_projected_vsd_region_local_rerank_loso_summary.json`

### 5.1 LOSO local rerank

| Route | Valid n | Exact Top1 | Exact Top3 | Median true rank |
|---|---:|---:|---:|---:|
| logCPM Top3 Network + local genes | 814 | 0.2101 | 0.4263 | 5.0 |
| native VSD Top3 Network + local genes | 814 | 0.2088 | 0.4238 | 5.0 |
| projected VSD Top3 Network + local genes | 814 | 0.2162 | 0.4337 | 4.0 |

正向结果：

- 加入 Network Top3 gate + local discriminative genes 后，exact Top3 从 direct exact 的 `0.3563` 提升到 `0.4337`。
- 说明 region-local rerank 是必要改进。

### 5.2 LOMO local rerank

| Route | Valid n | Exact Top1 | Exact Top3 | Median true rank |
|---|---:|---:|---:|---:|
| logCPM Top3 Network + local genes | 812 | 0.2057 | 0.3867 | 6.0 |
| native VSD Top3 Network + local genes | 812 | 0.2118 | 0.4027 | 5.0 |
| projected VSD Top3 Network + local genes | 812 | 0.2081 | 0.3978 | 5.0 |

限制：

- 在 LOMO local exact 层，projected VSD 并非最高；native VSD Top3 local route 的 exact Top3 更高。
- 这支持最终 hybrid：Network beam 使用 projected VSD，exact/local 层不要继续使用 pure projected VSD。

### 5.3 LOSO hybrid local rerank

| Route | Valid n | Exact Top1 | Exact Top3 | Median true rank |
|---|---:|---:|---:|---:|
| hybrid projected Network + logCPM local genes | 814 | 0.2088 | 0.4263 | 4.5 |
| logCPM Top3 Network + local genes | 814 | 0.2101 | 0.4263 | 5.0 |
| projected VSD Top3 Network + local genes | 814 | 0.2162 | 0.4337 | 4.0 |

解释：

- 单纯 local exact rerank 的 LOSO 结果并不足以单独证明 hybrid 最优。
- hybrid 的优势需要放在完整三级路线中解释：projected VSD 稳定生成 Network beam，logCPM downstream 提供更合理的 group/exact 层解释。

## 6. 完整正式三级内部验证

证据文件：

- LOSO: `results/bo2023_reference_projection_20260616_cleaned_symbols/formal_three_tier_loso_hybrid/hybrid_formal_loso_summary.json`
- LOMO: `results/bo2023_reference_projection_20260616_cleaned_symbols/formal_three_tier_lomo_hybrid/formal_lomo_validation_summary.json`

### 6.1 完整三级 LOSO hybrid

测试设计：complete three-tier LOSO；projected VSD Network Top3 beam；logCPM resolution/local exact rerank。

| 层级 | n | Top1 | Top3 | Median true rank |
|---|---:|---:|---:|---:|
| Network | 814 | 0.5835 | 0.9238 | 1.0 |
| Resolution group | 814 | 0.4447 | 0.7236 | 2.0 |
| Exact region | 814 | 0.2248 | 0.4533 | 4.0 |

额外事实：

- low-resolution predictions: `804/814 = 0.9877`
- high-resolution predictions: `10`
- conditional exact Top3 given Network Top1: `0.5516`

正向结果：

- Network Top3 `0.9238` 支持 Network beam 作为主 endpoint。
- Group Top3 `0.7236` 明显高于 Exact Top3 `0.4533`，支持使用 resolution group 作为 region-level 主报告粒度。
- Exact Top3 比早期 direct exact 和 local rerank baseline 更高，是完整三级路线的正向结果。

### 6.2 完整三级 LOMO hybrid

测试设计：formal three-tier leave-one-monkey-out；route-specific network、resolution-group 和 exact-region components fold-local 重建。

| 层级 | n | Top1 | Top3 | Median true rank |
|---|---:|---:|---:|---:|
| Network | 819 | 0.5775 | 0.9121 | 1.0 |
| Resolution group | 812 | 0.4138 | 0.6909 | 2.0 |
| Exact region | 812 | 0.2217 | 0.4236 | 5.0 |

额外事实：

- Network macro-by-monkey Top3: `0.8851`
- Group macro-by-monkey Top3: `0.6714`
- Exact macro-by-monkey Top3: `0.3911`
- low-resolution predictions: `809/812 = 0.9963`

正向结果：

- LOMO 中 Network Top3 仍为 `0.9121`，支持跨个体泛化。
- LOMO exact Top3 `0.4236` 与 LOSO `0.4533` 同一量级，说明 exact 层限制主要来自区域分辨率和相似脑区，而不仅仅是样本级泄漏。

不要写：

- “exact-region 已达到临床可直接判定。”
- “LOMO exact-region Top3 接近 Network Top3。”

## 7. AHBA 外部验证

证据文件：

- 早期 projected/logCPM 对比：`results/bo2023_reference_projection_20260616_cleaned_symbols/ahba_external_projected_vsd/ahba_projected_vsd_external_summary.json`
- 正式三级：`results/bo2023_reference_projection_20260616_cleaned_symbols/ahba_external_formal_three_tier/ahba_formal_three_tier_summary.json`
- special labels: `results/bo2023_reference_projection_20260616_cleaned_symbols/ahba_external_formal_three_tier/ahba_formal_three_tier_special_labels.csv`

### 7.1 早期 AHBA projected vs logCPM 对比

| Route | Network Top1 | Network Top3 | Exact Top1 | Exact Top3 |
|---|---:|---:|---:|---:|
| logCPM baseline | 0.5966 | 0.8283 | 0.2088 | 0.3297 |
| projected VSD | 0.4249 | 0.8884 | 0.0989 | 0.1648 |

解释：

- projected VSD 的 Network Top3 高于 logCPM baseline。
- projected VSD 的 exact-region 指标低于 logCPM baseline。
- 这正是后续 hybrid 路线的动机。

### 7.2 AHBA 正式三级验证

数据集：AHBA human RNA-seq raw counts。总样本 `242`；supported-for-accuracy 样本 `233`；exact-region evaluable 样本 `91`。

| Route | Network Top1 | Network Top3 | Group Top1 | Group Top3 | Exact Top1 | Exact Top3 |
|---|---:|---:|---:|---:|---:|---:|
| hybrid projected Network + logCPM exact | 0.7468 | 0.9442 | 0.3626 | 0.6703 | 0.2418 | 0.4286 |
| logCPM baseline | 0.6266 | 0.9700 | 0.2637 | 0.6264 | 0.1758 | 0.3077 |
| projected VSD only | 0.7468 | 0.9442 | 0.2637 | 0.5385 | 0.1099 | 0.2967 |

正向结果：

- hybrid 的 Network Top1 与 projected VSD 持平，高于 logCPM baseline。
- hybrid 的 Group Top3 和 Exact Top3 均高于 logCPM baseline 和 projected VSD only。
- 这支持“projected VSD 保住 Network 候选优势，logCPM downstream 修复 exact/group 层损失”的主线说法。

重要边界：

- AHBA 是 human normal brain RNA-seq，Bo2023 是 macaque atlas。
- exact 指标只对有稳定 Bo2023 映射的 AHBA 标签计算，不能泛化为所有人脑 exact-region。

### 7.3 AHBA 特定标签

hybrid route special labels：

| Label | n | Network Top3 | Group Top3 | Exact Top1 | Exact Top3 |
|---|---:|---:|---:|---:|---:|
| Insula | 8 | 1.0000 | 1.0000 | 0.3750 | 0.5000 |
| Caudate | 8 | 1.0000 | 1.0000 | 0.0000 | 0.6250 |
| Putamen | 9 | 1.0000 | 0.8889 | 0.3333 | 0.8889 |

写作建议：

- 可写 Putamen/Caudate/Insula 的 Network 层很稳。
- Putamen 的 exact Top3 较强；Caudate exact Top1 为 0，需要谨慎，只能强调 group/Top3。
- 不要把 Caudate 写成 exact Top1 成功案例。

## 8. TCGA/BraTS MRI-labeled 外部验证

证据文件：

- `results/bo2023_reference_projection_20260616_cleaned_symbols/tcga_labeled_hybrid_formal_external/tcga_labeled_hybrid_formal_summary.json`

数据集：TCGA/BraTS glioma RNA-seq with corrected MRI-derived labels。

- expression samples: `65`
- patients: `65`
- MRI truth exact labels are human atlas labels, not Bo2023 macaque region IDs。

| Route | Network Top1 | Network Top3 | Lobe Top1 | Lobe Top3 | Broad Top1 | Broad Top3 |
|---|---:|---:|---:|---:|---:|---:|
| hybrid projected Network + logCPM exact | 0.1538 | 0.4000 | 0.1077 | 0.2462 | 0.1385 | 0.6462 |
| logCPM baseline | 0.0923 | 0.3231 | 0.1385 | 0.2308 | 0.0615 | 0.2462 |
| projected VSD only | 0.1538 | 0.4000 | 0.2154 | 0.3692 | 0.1385 | 0.6462 |

正向结果：

- hybrid Network Top3 `0.4000` 高于 logCPM baseline `0.3231`。
- hybrid Broad anatomy Top3 `0.6462` 高于 logCPM baseline `0.2462`。

边界：

- 不能报告 Bo2023 exact-region accuracy。
- 不能把 TCGA/BraTS 写成 cfRNA 外部验证；它是 glioma tissue RNA-seq + MRI label。
- 该结果只能支持 coarse anatomical consistency。

## 9. 主线代码与 artifact 状态

本线程已将正式三级路线接入主线算法：

- Network projector route: `core/network_tracing.py`
- Three-tier region route: `core/bo2023_region_tracing.py`
- UI route text/trigger: `app/pages/tracing_page.py`
- sample expression columns preserved: `data_processor.py`
- projector artifact: `data/models/bo2023_reference_projector_linear_full.npz`

主线推理步骤：

```text
uploaded sample
  -> gene/QC/preprocess
  -> logCPM/logTPM-compatible query
  -> projector maps query to Bo2023-like VSD
  -> SaleemNetworks Top3 beam
  -> restrict Bo2023 region candidates to Network Top3
  -> logCPM local resolution-group rerank
  -> logCPM local exact-region rerank
  -> report Network, resolution group, exact candidate list
```

主线代码 smoke test 已跑过一个数据库样本 `19R348`：

- Network output scale: `projected_vsd`
- Network overlap: `199/200 = 0.995`
- Region route: `projected_vsd_network_top3_logcpm_resolution_local_exact`
- Region reference source: `raw_featurecounts_logcpm`
- Query source in that sample: `stored_log_tpm_fallback`

## 10. 论文初稿可写的正向结果

建议在初稿中只突出以下正向结果：

1. Projector 在 sample-level 上能稳定重建 VSD-like expression pattern：
   - global Pearson `0.9961`
   - median sample Pearson `0.9969`

2. projected VSD 是 Network candidate generation 的正向路线：
   - LOSO Network Top3 `0.9158`
   - LOMO Network Top3 `0.9133`
   - 完整三级 LOSO hybrid Network Top3 `0.9238`
   - 完整三级 LOMO hybrid Network Top3 `0.9121`

3. direct exact-region 不是最佳主线：
   - direct exact LOSO projected Top3 `0.3563`
   - direct exact LOMO projected Top3 `0.3116`
   - 两者都明显低于完整三级 hybrid 的 exact Top3。

4. local rerank 和 resolution group 是必要改进：
   - 完整三级 LOSO Group Top3 `0.7236`
   - 完整三级 LOMO Group Top3 `0.6909`
   - Group Top3 明显高于 Exact Top3，说明 group-level endpoint 更稳。

5. AHBA 外部验证支持 hybrid：
   - hybrid Network Top3 `0.9442`
   - hybrid Group Top3 `0.6703`
   - hybrid Exact Top3 `0.4286`
   - hybrid Exact Top3 高于 logCPM baseline `0.3077` 和 projected VSD only `0.2967`。

6. TCGA/BraTS 只支持粗粒度外部一致性：
   - hybrid Network Top3 `0.4000`
   - hybrid Broad anatomy Top3 `0.6462`

## 11. 论文初稿必须避免的说法

不要写：

- “整个 Bo2023 atlas 已重新投影成新图谱。”
- “projected VSD 在 exact-region 层全面优于 logCPM。”
- “TCGA/BraTS 验证了 Bo2023 exact-region accuracy。”
- “Caudate exact Top1 表现很好。”
- “所有 suspicious gene symbols 都是错误基因。”
- “ENSMFAG fallback 会使溯源无效。”
- “Network Top1 是足够稳定的唯一结论。”

可以写：

- “The VSD-projected query was used only for network-level candidate generation.”
- “Downstream resolution-group and exact-region reranking was performed in logCPM-compatible local expression space.”
- “Resolution groups are a more defensible region-level endpoint than exact Top1 calls.”
- “AHBA supports cross-species mapped-label validation, whereas TCGA/BraTS supports only coarse anatomical consistency.”

## 12. 建议给另一个线程的初稿改进任务

请另一个线程按以下任务撰写初稿：

1. 把论文 endpoint 改成三级层级，而不是 exact-region 单一 endpoint。
2. Methods 中明确说明：训练参考 atlas 是原始 Bo2023 VSD/logCPM 参考；只有 query 被 projector 映射到 VSD-like space 用于 Network beam。
3. Results 先写 Network 层内部 LOSO/LOMO，再写 full three-tier，再写 AHBA，再写 TCGA/BraTS。
4. AHBA 小节要强调 hybrid 的正向结果：保持 projected Network 优势，同时修复 projected-only exact 损失。
5. TCGA/BraTS 小节标题应使用 `coarse anatomical consistency`，不要使用 `exact region validation`。
6. Discussion 中明确 exact-region Top1 的限制，说明大量 low-resolution predictions 是 atlas 分辨率/相邻脑区相似性的真实问题。
7. Supplementary/Extended Data 放 gene-symbol audit、projector QC、GSE189919 projection feasibility 和 special-label 表。

## 13. 推荐图表和文件

已生成的论文指导图：

- `docs/figures/latest_three_tier_validation_summary.png`

已生成的图文报告：

- `output/pdf/reference_projection_test_report_zh.pdf`

可直接引用的图：

- `results/bo2023_reference_projection_20260616_cleaned_symbols/ahba_external_formal_three_tier/ahba_formal_three_tier_accuracy.png`
- `results/bo2023_reference_projection_20260616_cleaned_symbols/tcga_labeled_hybrid_formal_external/tcga_labeled_hybrid_formal_accuracy.png`
- `results/bo2023_reference_projection_20260616_cleaned_symbols/figures/projector_gene_qc_distributions.png`
- `results/bo2023_reference_projection_20260616_cleaned_symbols/figures/internal_network_accuracy.png`

## 14. 溯源标签补充：7m 与 VIP

本线程还确认了两个 Bo2023 region 的 Network 映射：

- `7m`: `Parietal, and Parieto-occipital region`
- `VIP`: `Parietal, and Parieto-occipital region`

二者在当前 resolution group 中属于：

```text
Parietal, and Parieto-occipital region::44563 + 5 + 7a + 7b + 7m + 7op + LIPv + VIP
```

这可用于解释 parietal/posterior parietal exact-region 混淆，不应把 `7m` 与 `VIP` 的互换解读为完全不同 Network 的错误。
