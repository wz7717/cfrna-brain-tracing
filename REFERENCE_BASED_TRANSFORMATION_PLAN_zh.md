# 训练集参数投影 / Reference-Based Transformation 实施计划

## 目的

本文档用于指导后续线程在本项目中实现和评估一条
`reference-based transformation` 分析分支。

核心目标是测试：

```text
外部样本或 count-derived 样本
能否先被投影到 Bo2023-like VSD 表达空间
再进入现有脑区溯源模型
```

这是一条探索性验证分支。除非通过严格内部验证，否则不能替换当前锁定的生产主线。

当前锁定生产主线是：

```text
Bo2023 VSD reference
-> fold-selected 200 genes
-> Pearson correlation
-> Top-3 pairwise rescue
```

## 核心思想

不要让每个外部数据集自己生成一个互不相干的 VST/VSD 尺度。

更合理的做法是先在 Bo2023 训练数据中学习一个转换规则：

```text
Bo2023 raw count-derived expression
-> Bo2023 VSD + batch-removed expression
```

然后把这个规则应用到留出样本或外部 count 样本：

```text
external raw count
-> count-derived expression，通常是 logCPM
-> Bo2023-trained projector
-> projected Bo2023-like VSD
-> existing tracing model
```

端到端要验证的主链路是：

```text
验证数据集 raw count 文件
-> 统一 gene ID / ortholog / gene symbol
-> 计算 logCPM
-> 套用只由 Bo2023 训练样本拟合出的 VSD 参数投影器
-> 生成 projected Bo2023-like VSD matrix
-> 用现有 Bo2023 VSD reference 溯源
-> 和现有 TPM/logCPM/rank 路线对照
```

注意：验证数据集 raw count 只能被 transform，不能参与 projector 参数拟合。
这里生成的是 `projected Bo2023-like VSD matrix`，不是外部数据集自己的 native VSD，
也不是严格复刻 DESeq2 的原生 VSD。

第一版最现实的实现不是完整复刻 DESeq2 的内部 VST 参数，而是做一个经验投影器：

```text
paired Bo2023 counts + Bo2023 VSD
-> fold-local empirical projector
```

## 现有输入文件

Bo2023 原始 count：

```text
bo2023 data/mfas5_819samples_28415genes_featurecounts_counts.txt
```

Bo2023 VSD + batch removed：

```text
bo2023 data/mfas5_819samples_23605genes_vsd4_rmbatch.xls
```

Bo2023 样本元数据：

```text
bo2023 data/Information of sequenced samples_update_full878_filter819.xlsx
```

在重新实现之前，应优先检查已有 VSD 重建 / frozen VST 相关产物：

```text
results/bo2023_vsd_reconstruction/bo2023_frozen_vst_reference.rds
results/bo2023_vsd_reconstruction/best_reconstructed_vsd.tsv.gz
scripts/reconstruct_bo2023_vsd.R
scripts/prepare_bo2023_vsd_metadata.py
```

本地已有 raw count 的外部数据集：

```text
data/external_validation/GSE106804/GSE106804_Gene_counts.txt.gz
data/external_validation/GSE189919/GSE189919_count.csv.gz
data/external_validation/GSE228512/GSE228512_hiseq_counts.txt.gz
data/external_validation/GSE228512/GSE228512_novaseq_counts.txt.gz
data/tcga_brain_tumor_expression/tcga_gbm_lgg_primary_tumor_unstranded_counts_sample_sum.tsv
data/ahba_human_rnaseq/raw_zips/H0351_2001_rnaseq.zip
data/ahba_human_rnaseq/raw_zips/H0351_2002_rnaseq.zip
```

AHBA zip 文件中包含 `RNAseqCounts.csv`，但当前 AHBA 验证脚本使用的是 `RNAseqTPM.csv`。

Ivy GAP 当前本地主要是 FPKM/TPM，不是严格 raw count：

```text
data/ivy_gap_anatomic_rnaseq/ivy_gap_anatomic_structure_fpkm_gene_symbol_matrix.tsv
data/ivy_gap_anatomic_rnaseq/ivy_gap_anatomic_structure_tpm_gene_symbol_matrix.tsv
```

除非重新获得 raw count，否则不要声称 Ivy GAP 做了严格 VST 投影。

## 非目标

不要把 projected values 直接称作真实 Bo2023 VSD，除非内部验证已经支持。

不要使用测试样本拟合任何转换参数。

不要使用外部队列标签、MRI 标签、肿瘤位置标签、诊断标签或结局标签来调 projector。

不要在验证通过前把 projected VSD 设为生产默认路线。

## 第一版推荐投影器

建议先做一个简单、可审计的按基因映射：

```text
输入：Bo2023 logCPM per gene
目标：Bo2023 VSD_batch_removed per gene
模型：每个基因一个 robust linear regression 或 ordinary linear regression
```

对每个基因 `g`：

```text
VSD_g = a_g * logCPM_g + b_g
```

logCPM 从 count 计算：

```text
CPM = count / sample_library_size * 1,000,000
logCPM = log1p(CPM)
```

需要设计 fallback 规则：

```text
如果某个基因训练集中非零样本太少：
    使用 z-score / median fallback，或排除该基因

如果预测值过于极端：
    clip 到 Bo2023 训练 VSD 分位数范围，例如 0.5%-99.5%
```

第一版代码应和生产推理路线隔离。

## 可比较的其他投影器

第一版 linear projector 跑通后，再比较这些变体：

1. per-gene robust linear projector
2. per-gene quantile mapping
3. per-gene z-score-to-reference mapping
4. rank-percentile projector
5. frozen DESeq2/VST replay，只有在现有 R 产物支持时才做

rank-based 输出应视为稳健 baseline，不应称作 VSD。

## 必须遵守的 Fold-Local 原则

Bo2023 LOSO：

```text
for each held-out sample:
    train_samples = all Bo2023 samples except held-out sample
    只用 train_samples 拟合 projector
    只用 train_samples 构建或选择 reference
    held-out sample count -> train-only projector -> projected VSD
    run tracing
```

Bo2023 LOMO：

```text
for each held-out monkey:
    train_samples = all Bo2023 samples from other monkeys
    只用 train_samples 拟合 projector
    只用 train_samples 构建或选择 reference
    held-out monkey samples -> train-only projector -> projected VSD
    run tracing
```

外部数据集：

```text
projector 只能在 Bo2023 training/reference samples 上 fit
external samples 只能 transform
然后 run tracing
```

如果使用外部样本整体分布估计任何参数，必须标记为：

```text
transductive sensitivity analysis
```

不能标记为 locked validation。

## 建议新增文件

先新增文件，不要直接修改锁定生产路线：

```text
scripts/build_bo2023_reference_projector.py
scripts/run_bo2023_projected_vsd_loso.py
scripts/run_bo2023_projected_vsd_lomo.py
scripts/apply_projected_vsd_to_external_counts.py
core/reference_projection.py
```

建议结果目录：

```text
results/bo2023_reference_projection_YYYYMMDD/
```

建议模型产物命名：

```text
bo2023_reference_projector_linear_fold_{fold_id}.npz
bo2023_reference_projector_linear_full.npz
```

建议结果文件：

```text
projector_gene_parameters.csv
projector_qc_summary.json
bo2023_projected_vsd_loso_detail.csv
bo2023_projected_vsd_loso_summary.json
bo2023_projected_vsd_lomo_detail.csv
bo2023_projected_vsd_lomo_summary.json
external_projected_vsd_<dataset>_detail.csv
external_projected_vsd_<dataset>_summary.json
method_note_reference_projection.md
```

## 阶段 1：数据审计

任务：

1. 读取 Bo2023 raw count 矩阵。
2. 读取 Bo2023 VSD/batch-removed 矩阵。
3. 读取 Bo2023 样本元数据。
4. 确认 count 和 VSD 的样本 ID 是否匹配。
5. 确认 gene ID / gene symbol。
6. 必要时合并重复基因。
7. 取 count genes 与 VSD genes 的交集。
8. 输出缺失样本和缺失基因报告。

产物：

```text
results/bo2023_reference_projection_YYYYMMDD/data_audit_summary.json
results/bo2023_reference_projection_YYYYMMDD/common_gene_panel.csv
```

最低可接受条件：

```text
n_common_samples 接近 819
n_common_genes 足以覆盖当前模型基因
locked model genes 要么存在，要么明确报告缺失
```

## 阶段 2：训练 Projector

任务：

1. Bo2023 raw counts -> logCPM。
2. 对齐 logCPM 与 VSD 矩阵的 genes 和 samples。
3. 在训练样本上拟合 per-gene projector。
4. 保存 slope、intercept、residual SD、training quantiles、QC flags。
5. 为不稳定基因增加 fallback 逻辑。

建议每个基因输出这些 QC 字段：

```text
gene_symbol
n_train_samples
n_nonzero_count_samples
logcpm_mean
logcpm_sd
vsd_mean
vsd_sd
slope
intercept
r2
spearman_r
residual_sd
fallback_reason
clip_low
clip_high
```

训练折内应检查：

```text
correlation(projected_vsd, native_vsd)
MAE(projected_vsd, native_vsd)
per-gene R2 distribution
per-sample correlation distribution
```

## 阶段 3：Bo2023 内部 LOSO 验证

至少比较这些路线：

```text
native_vsd:
    held-out native Bo2023 VSD sample -> fold-local Bo2023 VSD reference

projected_vsd:
    held-out Bo2023 count -> fold-local projector -> projected VSD -> fold-local Bo2023 VSD reference

logcpm_baseline:
    held-out Bo2023 logCPM -> compatible reference or direct correlation baseline

rank_baseline:
    held-out Bo2023 rank-percentile vector -> rank-transformed reference
```

主要指标：

```text
Network Top1
Network Top3
Exact Region Top1，如果支持
Exact Region Top3，如果支持
median true rank
abstain rate
n_overlap_genes
decision margin
```

判断标准：

```text
Projected VSD 应保留 native VSD 的大部分 Network Top3 性能。
Projected VSD 应不差于 plain logCPM / rank baseline。
如果 projected VSD 明显差于 logCPM，则不要继续做外部主张。
```

## 阶段 4：Bo2023 内部 LOMO 验证

使用 leave-one-monkey-out split 重复 projected VSD 验证。

LOMO 比 LOSO 更能说明泛化能力，因为留出的单位是一整只猴子，而不是单个样本。

必须记录：

```text
heldout_monkey_id
n_train_samples
n_test_samples
n_train_regions
n_test_regions
```

比较：

```text
native_vsd LOMO
projected_vsd LOMO
logCPM / rank baselines
```

## 阶段 5：外部 Count 数据集投影

只有阶段 3 通过后再做。

推荐顺序：

1. GSE189919 CSF RNA
2. GSE228512 EV RNA
3. GSE106804 tumor-specific EV RNA
4. AHBA RNA-seq counts
5. TCGA GBM/LGG counts

每个数据集需要：

1. 读取 raw counts。
2. 转换 gene ID 到 gene symbol；必要时做人-猕猴 ortholog 映射。
3. 合并重复 symbol。
4. count -> logCPM。
5. 应用 full Bo2023-trained projector。
6. 输出 projected Bo2023-like VSD matrix。
7. 运行现有 network tracing。
8. 和现有 TPM/logCPM/rank 路线对照。
9. 输出 confidence、margin、top-k distribution、overlap-gene QC。

除非有真实解剖标签，否则外部结果只能描述为：

```text
cross-domain stress test
```

不能描述为定位准确率验证。

## 阶段 6：报告

生成方法说明：

```text
results/bo2023_reference_projection_YYYYMMDD/method_note_reference_projection.md
```

必须写清：

```text
projector 使用 paired Bo2023 raw count-derived logCPM 和 Bo2023 VSD_batch_removed 训练。
所有内部验证均为 fold-local。
held-out sample 没有参与 projector 参数拟合。
外部 projected 结果不是 native Bo2023 VSD，应解释为 cross-domain projected-space analysis。
```

## 时间估计

最低可用版本：

```text
3-5 个工作日
```

范围：

```text
Bo2023 data audit
linear per-gene projector
Bo2023 LOSO validation
1 个外部数据集 proof-of-concept
```

研究级版本：

```text
7-10 个工作日
```

范围：

```text
LOSO + LOMO
多个 projector 变体
GSE106804、GSE189919、GSE228512、AHBA、TCGA
QC reports 和 method note
```

补充材料级版本：

```text
10-15 个工作日
```

范围：

```text
严格 fold-local 实现
ortholog audit
外部队列 QC
敏感性分析
论文级表格和图
```

## 主要风险

实现风险：

```text
Bo2023 count 和 VSD 的 sample IDs 可能需要仔细标准化。
count 矩阵基因数多于 VSD 矩阵。
VSD 矩阵已经 batch removed，简单 logCPM-to-VSD projector 未必能完整复现 batch removal。
人源外部数据需要 ortholog mapping。
肿瘤和液体活检表达分布与猕猴脑组织差异很大。
```

解释风险：

```text
Projected values 是近似 Bo2023-like VSD，不是真正 native VSD。
外部结果可能反映 domain shift、肿瘤生物学或生物流体成分，而不是解剖来源。
ComBat 或 transductive quantile mapping 如果不是 fold-local，会泄漏测试集信息。
```

## 后续线程建议起步命令

先检查这些文件：

```powershell
Get-ChildItem "bo2023 data"
Get-Content scripts\reconstruct_bo2023_vsd.R -TotalCount 220
Get-Content scripts\run_bo2023_loso_validation.py -TotalCount 260
Get-Content scripts\run_bo2023_leave_one_monkey_out_validation.py -TotalCount 260
Get-Content core\network_tracing.py -TotalCount 260
```

在构建任何模型前，先运行一个小型 data audit。

## 决策门槛

不要在以下条件满足前进入大规模外部验证：

```text
Bo2023 projected-VSD LOSO Network Top3 接近 native-VSD LOSO Network Top3；
并且 projected-VSD 不差于 simple logCPM / rank baseline。
```

如果这个门槛失败，应把 projected VSD 记录为 negative sensitivity analysis，而不是模型改进。
