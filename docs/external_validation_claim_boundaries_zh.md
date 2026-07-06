# 外部验证数据集的真值边界与可声明结论

本文档专门解释 cfRNA-BrainTrace v0.1.6 的外部验证为什么分成三类：

```text
AHBA       = mapped-label external validation
TCGA/BraTS = coarse anatomical consistency
GSE189919  = projection feasibility only
```

如果只记一句话：

```text
外部验证能声明什么，取决于外部数据有没有解剖真值、真值分辨率有多细、真值标签能否合理映射到 Bo2023 macaque-derived hierarchy；不能因为模型能输出预测，就把没有真值的数据当成 accuracy validation。
```

## 1. 为什么外部验证要先讲 truth label

一个验证结果能不能叫 accuracy，取决于是否有可比较的 truth label。

例如：

```text
预测: Network A
真值: Network A
=> 可以算 hit
```

但如果没有真值：

```text
预测: Network A
真值: unknown
=> 不能算对，也不能算错
```

因此外部验证第一步不是跑模型，而是问：

```text
这个数据集有什么 truth?
truth 的分辨率是什么?
truth 能否映射到当前 hierarchy?
哪些样本支持 accuracy?
哪些样本只能做 feasibility 或 stress test?
```

## 2. 当前外部验证总览

当前 Table S4 包含三类外部验证：

| Dataset | Truth label type | 可以声明什么 | 不可以声明什么 |
|---|---|---|---|
| AHBA | human brain anatomical labels mapped to hierarchy | mapped-label Network/group/exact support | 直接 human-to-macaque exact anatomical identity |
| TCGA/BraTS | human MRI-derived tumour labels | coarse anatomical consistency | Bo2023 macaque exact-region localization |
| GSE189919 | no patient-level anatomical truth | projection feasibility / gene-space coverage | localization accuracy |

对应 Table S4 数字：

```text
AHBA Network Top1/Top3 = 74.68% / 94.42%
AHBA Resolution group Top1/Top3 = 36.26% / 67.03%
AHBA Exact region Top1/Top3 = 24.18% / 42.86%

TCGA/BraTS Network Top1/Top3 = 15.38% / 40.00%
TCGA/BraTS Broad anatomy Top1/Top3 = 13.85% / 64.62%

GSE189919 projector gene overlap = 15622 / 21668
GSE189919 coverage = 72.10%
```

## 3. AHBA 是什么验证

### 3.1 数据集性质

AHBA 是 Allen Human Brain Atlas human brain RNA-seq。

它的优点：

- 是正常 human brain expression。
- 有 anatomical sample labels。
- 可以把部分 human anatomical labels 映射到当前 Network / group / exact hierarchy。

它的限制：

- 它是 human，不是 macaque。
- AHBA label 和 Bo2023 region label 不是天然一一对应。
- 有些 label 只能映射到粗层级，有些能映射到 exact-evaluable subset。

所以 AHBA 是：

```text
mapped-label external validation
```

不是：

```text
直接同物种 exact-region truth validation
```

### 3.2 主要脚本

```text
scripts/run_ahba_projected_vsd_formal_three_tier_external.py
```

相关辅助：

```text
scripts/run_ahba_human_rnaseq_external_validation.py
scripts/run_ahba_projected_vsd_external_validation.py
scripts/export_ahba_external_validation_figures.py
```

### 3.3 脚本怎样处理 AHBA truth

核心逻辑是：

1. 读取 AHBA expression 和 metadata。
2. 读取 AHBA 的 `main_structure`、`sub_structure` 等 anatomical labels。
3. 调用 label mapping 逻辑：

```text
mapping_for_ahba_label(...)
```

4. 为每个 AHBA sample 标记：

```text
supported_for_accuracy
accuracy_level
allowed_networks
allowed_regions
```

5. 只在支持的层级上计算 accuracy。

这意味着：

- 如果只支持 Network mapping，就只能算 Network。
- 如果支持 exact mapping，才进入 exact-evaluable subset。

### 3.4 AHBA 的 route 设置

脚本中比较多个 route，其中当前投稿重点是 hybrid route：

```text
projected-VSD Network beam
logCPM-compatible group/exact reranking
```

还报告两个对照：

```text
logCPM baseline exact
projected-VSD-only exact
```

对照的目的不是替代主路线，而是说明 hybrid route 在 exact Top3 上优于这些 side routes。

### 3.5 AHBA 输出文件

主要输出：

```text
ahba_formal_three_tier_sample_detail.csv
ahba_formal_three_tier_metrics.csv
ahba_formal_three_tier_special_labels.csv
ahba_formal_three_tier_resolution_audit.csv
ahba_formal_three_tier_summary.json
ahba_formal_three_tier_summary.md
ahba_formal_three_tier_accuracy.png
```

### 3.6 AHBA 当前论文数字

当前 Table S4：

```text
Hybrid Network:
  Top1 = 74.68%
  Top3 = 94.42%

Hybrid Resolution group:
  Top1 = 36.26%
  Top3 = 67.03%

Hybrid Exact region:
  Top1 = 24.18%
  Top3 = 42.86%

logCPM baseline Exact:
  Top1 = 17.58%
  Top3 = 30.77%

Projected-VSD-only Exact:
  Top1 = 10.99%
  Top3 = 29.67%
```

Supplement 中还说明 AHBA 样本范围：

```text
242 total
233 supported
91 exact-evaluable
```

### 3.7 AHBA 能说什么

可以说：

```text
AHBA supports cross-species mapped-label transfer.
Hybrid route shows strong Network mapped-label support.
Exact-region result is only for exact-evaluable mapped labels.
```

中文表述：

```text
AHBA 支持外部 mapped-label 转移验证，尤其是 Network 层较稳健；exact-region 只限于能稳定映射的 AHBA 标签子集。
```

### 3.8 AHBA 不能说什么

不能说：

```text
AHBA 证明模型可以在人脑中精确定位 macaque exact region。
```

原因：

- AHBA 是 human brain。
- Bo2023 hierarchy 来自 macaque reference。
- human label 到 macaque region 是 harmonized mapping，不是直接解剖同一性。

## 4. TCGA/BraTS 是什么验证

### 4.1 数据集性质

TCGA/BraTS 这里使用的是 glioma tissue RNA-seq 和 MRI-derived labels。

它的优点：

- 有 human MRI-derived spatial/coarse labels。
- 可以测试模型输出在肿瘤组织上的粗解剖一致性。

它的限制非常重要：

- 样本是 glioma tumour tissue，不是正常脑组织。
- truth label 来自 human MRI/tumour context。
- truth label 不是 Bo2023 macaque exact-region identifier。

所以 TCGA/BraTS 是：

```text
coarse anatomical consistency
```

不是：

```text
macaque exact-region validation
```

### 4.2 主要脚本

```text
scripts/run_tcga_labeled_hybrid_formal_external.py
```

相关脚本：

```text
scripts/build_corrected_brats_mri_truth.py
scripts/evaluate_brats_tcga_lgg_65_mri_truth.py
scripts/analyze_tcga_gbm_lgg_network_domain_shift.py
```

### 4.3 脚本怎样处理 TCGA/BraTS truth

脚本读取：

```text
--tcga-counts
--manifest
--labels
--bo-counts
--bo-vsd
--sample-info
--region-meta
```

核心函数：

```text
split_candidates
norm
hit
candidate_regions
summarize
```

其中 `hit(predicted, truth)` 做的是：

```text
预测候选列表与 truth candidates 是否有交集
```

这适合粗标签一致性，因为 MRI-derived truth 可能不是单一精确 Bo2023 region。

### 4.4 TCGA/BraTS 输出层级

脚本会输出多个层级：

```text
Network
Lobe
Broad anatomy
```

但论文当前 Table S4 采用：

```text
Network
Broad anatomy
```

并明确写：

```text
Coarse consistency only
```

### 4.5 TCGA/BraTS 输出文件

主要输出：

```text
tcga_labeled_hybrid_formal_sample_detail.csv
tcga_labeled_hybrid_formal_metrics.csv
tcga_labeled_hybrid_formal_summary.json
tcga_labeled_hybrid_formal_summary.md
tcga_labeled_hybrid_formal_metrics.png
```

### 4.6 TCGA/BraTS 当前论文数字

当前 Table S4 route：

```text
n = 65

Hybrid Network:
  Top1 = 15.38%
  Top3 = 40.00%

Hybrid Broad anatomy:
  Top1 = 13.85%
  Top3 = 64.62%
```

### 4.7 TCGA/BraTS 能说什么

可以说：

```text
TCGA/BraTS supports coarse anatomical consistency in human glioma tissue.
Broad anatomy Top3 is 64.62%.
```

中文表述：

```text
TCGA/BraTS 只支持肿瘤组织中的粗解剖一致性，尤其是 broad anatomy 层面的候选一致性。
```

### 4.8 TCGA/BraTS 不能说什么

不能说：

```text
TCGA/BraTS 验证了 macaque Network-level localization。
TCGA/BraTS 验证了 Bo2023 exact-region localization。
TCGA/BraTS 可以直接证明临床定位能力。
```

原因：

- MRI truth 是 human imaging/tumour labels。
- tumour expression 有 disease/domain shift。
- truth label 粒度与 Bo2023 macaque exact region 不等价。

### 4.9 旧 MRI side route 数字

曾出现另一条 MRI truth / matching rule side route：

```text
Network Top3 = 36.92%
Broad Top3 = 80.00%
```

这不是当前 Table S4 submission route。

正确处理：

```text
legacy / inconsistency / alternate matching rule
```

不能替换当前投稿数字：

```text
Network Top3 = 40.00%
Broad anatomy Top3 = 64.62%
```

## 5. GSE189919 是什么验证

### 5.1 数据集性质

GSE189919 是外部 CSF / biofluid expression matrix。

关键限制：

```text
没有 patient-level anatomical truth
```

这意味着它不能判断：

```text
预测 Network 是否正确
预测 region 是否正确
```

因此它不是 accuracy validation。

### 5.2 主要脚本

```text
scripts/analyze_gse189919_csf_tracing_validation.py
```

### 5.3 脚本读取什么

默认读取：

```text
GSE189919_count.csv.gz
GSE189919_tpm_count.csv.gz
GSE189919_family.soft.gz
```

入口参数：

```text
--data-dir
--outdir
```

### 5.4 GSE189919 脚本做什么

它可以做：

1. 检查 gene overlap。
2. 检查输入矩阵是否可投影。
3. 生成 Network prediction distributions。
4. 计算不同 route 的 agreement。
5. 输出 sample QC。
6. 做 methodological review。

它不能做：

```text
accuracy
Top1 hit
Top3 hit against anatomical truth
```

因为没有 truth。

### 5.5 GSE189919 输出文件

常见输出：

```text
gse189919_sample_predictions.csv
sample_input_qc.csv
sample_detection_summary.csv
network_broad_top1_occurrence.csv
network_broad_top3_occurrence.csv
network_broad_prediction_collapse.csv
route_agreement.csv
algorithm_audit.csv
methodological_review.csv
```

### 5.6 GSE189919 当前论文数字

当前 Table S4：

```text
Projector gene overlap = 15622 / 21668
Coverage = 72.10%
Conclusion = No accuracy claim
```

可以写成：

```text
GSE189919 supports projection feasibility because 15,622 of 21,668 projector genes were present.
```

不能写成：

```text
GSE189919 validated brain-source localization accuracy.
```

### 5.7 algorithm audit failed 怎么解释

GSE189919 脚本里有 optional algorithm audit。

如果看到：

```text
Algorithm audit failed:
baseline_production_tpm ...
baseline_production_count_cpm ...
```

这不表示当前投稿准确率路线失败。

正确解释是：

```text
旧 baseline / production comparison 的某些技术假设不应作为当前 submission accuracy evidence。
```

workflow 中应标为：

```text
EXPECTED_LEGACY_ALGORITHM_LIMITATION
```

而不是把 GSE189919 当定位准确率路线。

## 6. 三个外部数据集的 allowed claims

### 6.1 AHBA allowed claims

允许：

```text
mapped-label external validation
cross-species mapped-label transfer
Network-level support
group-level support
exact metrics only for exact-evaluable mapped labels
```

不允许：

```text
direct human-macaque exact anatomical equivalence
clinical localization proof
unmapped labels as accuracy evidence
```

### 6.2 TCGA/BraTS allowed claims

允许：

```text
coarse anatomical consistency
tumour-tissue stress test
broad anatomy candidate consistency
```

不允许：

```text
Bo2023 exact-region validation
normal brain source localization accuracy
clinical diagnostic localization
```

### 6.3 GSE189919 allowed claims

允许：

```text
projection feasibility
gene-space coverage
biofluid transfer stress test
prediction distribution audit
```

不允许：

```text
accuracy
Top1/Top3 localization hit
brain-source truth validation
```

## 7. 为什么不能“模型有输出就算验证”

模型对任何输入都可能输出一个 ranking。

但 ranking 本身不等于验证。

验证需要：

```text
prediction
truth
comparison rule
denominator
claim boundary
```

如果缺少 truth，只能做：

```text
feasibility
QC
distribution audit
stress test
```

GSE189919 就属于这一类。

## 8. 怎样读 Table S4

Table S4 的每一行都要同时读：

```text
Dataset
Route
Endpoint
Top1
Top3
Conclusion
```

不能只摘 Top3 数字。

例如：

```text
TCGA/BraTS Broad anatomy Top3 = 64.62%
Conclusion = Coarse consistency only
```

如果只写：

```text
TCGA/BraTS Top3 = 64.62%
```

就会误导读者以为这是 exact localization accuracy。

## 9. 审稿复现时的检查点

### 9.1 AHBA

检查：

```text
ahba_formal_three_tier_metrics.csv
ahba_formal_three_tier_sample_detail.csv
ahba_formal_three_tier_special_labels.csv
```

确认：

```text
supported_for_accuracy
accuracy_level
exact-evaluable subset
```

数字应对应：

```text
Network Top3 = 94.42%
Resolution group Top3 = 67.03%
Exact region Top3 = 42.86%
```

### 9.2 TCGA/BraTS

检查：

```text
tcga_labeled_hybrid_formal_metrics.csv
tcga_labeled_hybrid_formal_summary.md
```

确认 summary 中有类似边界说明：

```text
Exact Bo2023 region accuracy is not reported because MRI truth regions use human atlas labels.
```

数字应对应：

```text
Network Top3 = 40.00%
Broad anatomy Top3 = 64.62%
```

### 9.3 GSE189919

检查：

```text
sample_input_qc.csv
route_agreement.csv
algorithm_audit.csv
methodological_review.csv
```

确认最终解释是：

```text
projection feasibility only
no localization accuracy calculation
```

数字应对应：

```text
15622 / 21668
72.10%
```

## 10. 最短复述

外部验证的边界可以这样复述：

```text
AHBA 有 human anatomical labels，可映射到当前 hierarchy，所以能做 mapped-label validation，但不是直接 human-macaque exact identity。
TCGA/BraTS 有 MRI-derived tumour labels，只能支持 coarse anatomical consistency，不能验证 Bo2023 exact regions。
GSE189919 没有 patient-level anatomical truth，只能说明 gene overlap 和 projection feasibility，不能算 accuracy。
```
