# Network Beam 算法深入讲解

本文档专门解释 cfRNA-BrainTrace v0.1.6 的第一层算法：Network Top3 beam 是怎样从表达矩阵算出来的。

如果只记一句话：

```text
用户表达矩阵先变成 logCPM，再用 projector 映射到 Bo2023-like VSD 空间，然后与每个 Saleem Network 的 reference centroid 做 Pearson correlation，取分数最高的前三个 Network 作为 Top3 beam。
```

这个 Top3 beam 是后续 resolution-group 和 exact-region reranking 的候选范围。projected-VSD 的主要用途就在这里。

## 1. Network beam 在整个项目里的位置

完整投稿路线是：

```text
raw counts / logCPM query
        |
        v
logCPM-compatible query vector
        |
        v
logCPM-to-VSD projector
        |
        v
projected-VSD query vector
        |
        v
Network centroid correlation
        |
        v
Network Top3 beam
        |
        v
logCPM-compatible local reranking
        |
        +--> resolution-group candidates
        +--> exact-region exploratory candidates
```

因此 Network beam 是一个“先缩小范围”的步骤。

它不直接给最终 exact brain region，而是回答：

> 这个样本最像哪几个粗粒度脑网络？

论文里 Network Top3 是最稳健的 endpoint。region 层只在这个 beam 内继续排序。

## 2. 相关代码文件

Network beam 相关代码主要在这些文件：

```text
core/reference_projection.py
core/network_tracing.py
scripts/run_bo2023_projected_vsd_loso.py
scripts/run_bo2023_projected_vsd_lomo.py
scripts/run_bo2023_hybrid_formal_loso.py
scripts/run_bo2023_projected_vsd_formal_lomo.py
```

其中：

- `core/reference_projection.py`
  - 负责 raw counts 到 logCPM、logCPM 到 projected-VSD。

- `core/network_tracing.py`
  - 负责 app/生产环境的 Network ranking。

- `run_bo2023_projected_vsd_loso.py`
  - 验证 projected-VSD Network beam 的 LOSO 表现。

- `run_bo2023_projected_vsd_lomo.py`
  - 验证 projected-VSD Network beam 的 LOMO 表现。

- `run_bo2023_hybrid_formal_loso.py`
  - 完整 three-tier LOSO 中的 Network 层。

- `run_bo2023_projected_vsd_formal_lomo.py`
  - 完整 three-tier LOMO 中的 Network 层。

## 3. 输入数据长什么样

### 3.1 用户或外部样本输入

最理想输入是 gene-level raw counts：

```text
gene_symbol  read_count
GAD1         123
SLC17A7      456
MBP          789
```

也可以是已经处理过的 logCPM-like 表达。

TPM/logTPM 在 app 中也可接受，但它是 legacy fallback，不应等同于当前 Bo2023 logCPM 验证路线。

### 3.2 Bo2023 reference

Bo2023 reference 有三类关键数据：

1. raw counts
   - 用于计算 logCPM。

2. native VSD matrix
   - Bo2023 内部已有的 variance-stabilized expression。

3. sample metadata
   - 每个样本对应的 `Region`、`SaleemNetworks`、`MonkeyID`。

### 3.3 Network model

生产 Network model 文件：

```text
data/models/bo2023_saleem_network_top200_model.npz
```

里面主要包含：

```text
genes
networks
reference
fisher_scores
```

可以把 `reference` 理解为：

```text
每个 Network 的代表性表达向量
```

即 Network centroid。

## 4. 第一步：raw counts 变成 logCPM

代码位置：

```text
core/reference_projection.py
```

函数：

```python
compute_logcpm(counts)
```

公式：

```text
library_size = 一个样本所有 gene counts 的总和
CPM_gene = count_gene / library_size * 1,000,000
logCPM_gene = log1p(CPM_gene)
```

为什么要这么做：

- raw counts 受测序深度影响很大。
- CPM 先把不同样本归一到每百万 reads。
- log1p 压缩高表达基因的极端值。

直观例子：

```text
样本 A 总 reads = 10,000,000
某 gene count = 1,000
CPM = 1,000 / 10,000,000 * 1,000,000 = 100
logCPM = log(1 + 100)
```

这样不同测序深度的样本更可比较。

## 5. 第二步：logCPM 映射到 projected-VSD

代码位置：

```text
core/reference_projection.py
```

核心函数：

```python
fit_linear_projector(...)
apply_projector(...)
```

### 5.1 为什么需要 projector

Bo2023 的 native VSD 对 Network beam 更稳定，但外部用户通常只有 raw counts 或 logCPM，没有 DESeq2/VSD 处理上下文。

如果要求用户上传 VSD，会降低可复现性，因为 VSD 通常依赖整个批次和 reference 处理方式。

所以项目训练了一个简单、可审计的 gene-wise projector：

```text
logCPM -> Bo2023-like VSD
```

### 5.2 projector 的数学形式

对每一个 gene 单独拟合：

```text
VSD_gene = slope_gene * logCPM_gene + intercept_gene
```

这不是复杂黑箱模型。每个 gene 只有两个主要参数：

```text
slope
intercept
```

此外还有：

```text
clip_low
clip_high
fallback_reason
```

### 5.3 为什么要 clip

外部 query 的表达值可能超出 Bo2023 训练分布。

如果直接线性外推，某些 gene 的 projected-VSD 可能离谱。代码会把 projected 值裁剪到训练集中 native VSD 的合理范围：

```python
projected = np.clip(projected, fit.clip_low[:, None], fit.clip_high[:, None])
```

直观理解：

> 如果训练集中某 gene 的 VSD 合理范围大概是 2 到 12，外部样本投影出 50，就把它拉回合理范围内。

### 5.4 fallback gene 是什么

某些 gene 在训练集中表达太低或 logCPM 方差太小，无法稳定拟合 slope。

这些 gene 不强行建线性模型，而是 fallback 到训练集 VSD 均值。

这样做牺牲了一点信息，但避免低质量 gene 引入噪声。

### 5.5 projector 文件

生产 projector 存在：

```text
data/models/bo2023_reference_projector_linear_full.npz
```

它包含：

```text
genes
slope
intercept
clip_low
clip_high
logcpm_mean
logcpm_sd
vsd_mean
vsd_sd
fallback_reason
```

## 6. 第三步：构建 Network centroid

Network centroid 是每个 Network 的代表性表达向量。

在验证脚本里，centroid 通常由训练折样本计算：

```text
Network A centroid = training samples in Network A 的平均表达
Network B centroid = training samples in Network B 的平均表达
...
```

代码中常见函数名：

```python
build_centroids(...)
build_label_reference(...)
```

### 6.1 为什么验证中要 fold-local centroid

LOSO/LOMO 验证时，不能把测试样本的信息泄漏到 reference 中。

因此每个 fold 都要重新用训练样本构建 centroid：

```text
LOSO:
  train = 818 samples
  test = 1 sample

LOMO:
  train = other monkeys
  test = held-out monkey
```

这叫 fold-local reference。

### 6.2 生产 app 与验证脚本的区别

生产 app 使用锁定好的 model：

```text
data/models/bo2023_saleem_network_top200_model.npz
```

验证脚本为了防止数据泄漏，会在每个 fold 重新构建训练 reference。

这两个场景目的不同：

- app：给用户一个锁定模型做预测。
- validation：证明模型路线在留出数据上能泛化。

## 7. 第四步：Pearson correlation 打分

Network scoring 的核心是 Pearson correlation。

对每个 Network：

```text
score(Network k) = corr(query_projected_VSD, centroid_Network_k)
```

相关代码：

```python
scores = trace_corr(model["reference"], vector)
```

或验证脚本里的：

```python
corr_scores(reference, sample)
correlation_scores(reference, sample)
```

### 7.1 为什么用 correlation

Pearson correlation 看的是表达模式是否相似，而不是绝对表达量是否一样。

这对跨样本、跨批次更稳健：

- 如果两个样本整体测序深度不同，但高低表达模式相似，correlation 仍然高。
- 对候选排序来说，模式相似性比绝对数值更重要。

### 7.2 correlation 怎么变成排名

假设有 10 个 Saleem Networks，得到 10 个分数：

```text
Network A: 0.71
Network B: 0.34
Network C: 0.55
...
```

按分数从高到低排序：

```text
Top1 = Network A
Top2 = Network C
Top3 = Network B
```

这三个就是 Network Top3 beam。

## 8. 第五步：Top1、Top3、true rank 怎么算

在验证脚本里，每个样本会生成一行：

```text
truth_network
pred_top1
pred_top2
pred_top3
hit1
hit3
true_rank
```

### 8.1 Top1 hit

```text
hit1 = pred_top1 == truth_network
```

### 8.2 Top3 hit

```text
hit3 = truth_network in [pred_top1, pred_top2, pred_top3]
```

### 8.3 true rank

如果真实 Network 在完整排序里第 2 名：

```text
true_rank = 2
```

如果真实 Network 不在 Top3，但在第 5 名：

```text
true_rank = 5
```

### 8.4 accuracy

所有样本汇总：

```text
Top1 accuracy = mean(hit1)
Top3 accuracy = mean(hit3)
```

所以分母就是 detail 表的行数。

## 9. pairwise rescue 是什么

代码位置：

```text
core/network_tracing.py
scripts/run_bo2023_projected_vsd_formal_lomo.py
```

相关函数：

```python
_apply_pairwise_rescue(...)
evaluate_pairwise_rescue(...)
```

### 9.1 为什么需要 rescue

有些 Network 在训练中经常互相混淆。全局 centroid correlation 有时会把两个非常接近的 Network 排反。

pairwise rescue 的思路是：

> 如果 Top1 和 Top2/Top3 是一对已知容易混淆的 Network，就用专门区分这两个 Network 的 pair-specific genes 再比较一次。

### 9.2 rescue 不做什么

它不重新搜索所有 Networks。

它不把 Top3 外的 Network 加进来。

它不改变“Top3 beam 是候选范围”的基本定义。

它只可能在原 Top3 内部调整 Top1：

```text
原排序: A, B, C
pairwise rescue 后: B, A, C
```

Top3 集合仍然是：

```text
{A, B, C}
```

### 9.3 为什么这点重要

因为后续 region reranking 依赖 Network Top3 beam 选 candidate regions。

如果 rescue 可以引入 Top3 外的新 Network，region 候选范围会变得不稳定。当前实现避免了这个问题。

## 10. Network-only validation 与 formal Network 的区别

项目里有两类 Network 数字：

1. projected-VSD Network-only validation
2. formal three-tier route 里的 Network layer

它们都看 Network，但目的不同。

### 10.1 projected-VSD Network LOSO/LOMO

脚本：

```text
scripts/run_bo2023_projected_vsd_loso.py
scripts/run_bo2023_projected_vsd_lomo.py
```

问题：

> projected-VSD 作为 Network beam 表示是否有效？

论文数字：

```text
LOSO:
  n = 819
  Top1 = 58.00%
  Top3 = 91.58%

LOMO:
  n = 819
  Top1 = 53.72%
  Top3 = 91.33%
```

这证明 projected-VSD 可以作为 Network beam 的基础。

### 10.2 formal LOSO/LOMO Network

脚本：

```text
scripts/run_bo2023_hybrid_formal_loso.py
scripts/run_bo2023_projected_vsd_formal_lomo.py
```

问题：

> 完整投稿路线的第一层 Network 表现如何？

论文数字：

```text
Formal LOSO Network:
  n = 819
  Top1 = 58.24%
  Top3 = 92.19%

Formal LOMO Network:
  n = 819
  Top1 = 57.75%
  Top3 = 91.21%
```

这些数字属于完整 three-tier route，因此要和后续 group/exact 的结果放在同一套 validation framework 里理解。

### 10.3 为什么不能混写

因为：

- Network-only validation 是方法合理性验证。
- formal Network 是投稿主流程中的第一层结果。
- formal LOMO 可能带有 pairwise Top1 rescue 等 formal route 设定。

所以：

```text
projected-VSD Network LOMO Top3 = 91.33%
formal LOMO Network Top3 = 91.21%
```

两者接近，但不是同一条路线。

## 11. formal LOSO 中 Network beam 怎么跑

代码位置：

```text
scripts/run_bo2023_hybrid_formal_loso.py
```

关键代码逻辑：

```text
network_reference = build_centroids(vsd_values[locked_rows, :], network_labels, networks, train_idx)
projected_locked = loo_project_rows(logcpm_values, vsd_values, locked_rows, sample_idx)
network_scores = correlation_scores(network_reference, projected_locked)
network_top = sort(network_scores)[0:3]
network_rows.append(...)
```

白话解释：

1. 用训练样本的 VSD 表达建立 Network reference。
2. 把 held-out sample 的 logCPM 投影到 VSD-like 空间。
3. 只用 locked Network genes。
4. 与每个 Network reference 做 correlation。
5. 排名前三就是 Network Top3。
6. 不管 region 是否可评价，Network 结果都先写入。

这就是为什么 formal LOSO Network 分母是 819。

## 12. formal LOMO 中 Network beam 怎么跑

代码位置：

```text
scripts/run_bo2023_projected_vsd_formal_lomo.py
```

关键逻辑：

```text
test_indices = held-out monkey samples
train_indices = other monkeys
projected_test = fit_project_rows(logcpm_values, vsd_values, all_genes, train_indices, test_indices)
network_reference = build_label_reference(...)
scores = correlation_scores(...)
pair_detail, pair_prob = evaluate_pairwise_rescue(...)
```

白话解释：

1. 每次留出整只 monkey。
2. 用其他 monkeys 训练 projection/reference。
3. 把 held-out monkey 的样本投影为 projected-VSD。
4. 与训练 monkeys 的 Network reference 做 correlation。
5. 得到 Top3。
6. 在 Top3 内做 pairwise rescue。
7. 写入 Network detail。

formal LOMO 的 Network 分母仍是 819，因为每个样本都有 Network truth label。

## 13. 生产 app 中 Network beam 怎么跑

代码位置：

```text
core/network_tracing.py
```

函数：

```python
trace_network_expression(...)
```

生产流程：

1. 用户上传 expression。
2. `_sample_logcpm_series` 判断输入尺度。
3. 如果启用 projector，就用 `apply_projector`。
4. 与锁定 Network model 的 reference 做 correlation。
5. 可选 pairwise rescue。
6. 输出排序结果和 metadata。

返回结构大致是：

```text
results:
  network_id
  rank
  score
  confidence

meta:
  endpoint
  method
  n_networks
  n_model_genes
  n_overlap_genes
  overlap_fraction
  traceability
  pairwise_rescue
  reference_projection
```

审稿或排错时要重点看：

```text
meta.reference_projection.enabled
meta.reference_projection.output_scale
meta.overlap_fraction
meta.pairwise_rescue
```

## 14. 为什么 projected-VSD 只用于 Network beam

内部验证显示：

- projected-VSD 对粗 Network Top3 很稳健。
- direct exact-region projected-VSD scoring 更不稳定，尤其在 LOMO。

所以当前路线把 projected-VSD 限定在第一层：

```text
projected-VSD -> Network Top3 beam
```

后续 region 层使用：

```text
logCPM-compatible local reranking
```

这避免了把 projected-VSD 过度解释成 exact-region 定位能力。

## 15. 常见问题

### 15.1 Network Top3 是不是最终答案

不是。它是第一层候选束。

最终 app 可以显示 Network、resolution-group 和 exact-region candidates，但论文主张强度不同：

- Network Top3 最稳健。
- Resolution group 是较可辩护的 region-level 输出。
- Exact region 是 exploratory ranking。

### 15.2 如果真实 Network 在 Top3，但不是 Top1，算成功吗

对 Top3 指标来说算成功，对 Top1 指标不算。

例如：

```text
truth = Network C
prediction = [Network A, Network B, Network C]
Top1 = fail
Top3 = hit
```

### 15.3 为什么 Top3 比 Top1 更重要

因为 cfRNA-BrainTrace 是候选排序工具，不是单标签诊断器。

Top3 更符合实际用途：

```text
给出最可能的几个来源范围，让后续 region reranking 和人工解释继续收窄。
```

### 15.4 pairwise rescue 会不会造成数据泄漏

在验证脚本中，pairwise rescue model 是 fold-local 构建或按训练设定约束的。它不能使用 held-out truth 来调整该样本结果。

并且 rescue 只在原 Top3 内部调整顺序，不改变 Top3 候选集合。

### 15.5 为什么 formal LOSO Network 没有因为 5 个 unsupported region 样本降到 814

因为 unsupported 是 region-level support 问题，不是 Network truth 缺失。

那 5 个样本仍有 Network label，也能合法评价 Network beam。

所以：

```text
Network n = 819
group/exact n = 814
```

## 16. 审稿复现时怎样确认 Network beam 路线没跑错

检查点：

1. projected-VSD Network LOSO 输出：

```text
bo2023_projected_vsd_loso_route_summary.csv
Top1 = 58.00%
Top3 = 91.58%
n = 819
```

2. projected-VSD Network LOMO 输出：

```text
bo2023_projected_vsd_lomo_route_summary.csv
Top1 = 53.72%
Top3 = 91.33%
n = 819
```

3. formal LOSO Network 输出：

```text
hybrid_formal_loso_network_route_metrics.csv
Top1 = 58.24%
Top3 = 92.19%
n = 819
```

4. formal LOMO Network 输出：

```text
formal_lomo_network_route_metrics.csv
Top1 = 57.75%
Top3 = 91.21%
n = 819
```

5. 论文和 workflow 中不能把 `92.38%` 当当前 Network Top3。

`92.38%` 是旧的 conditional denominator 数字，只能作为 legacy denominator inconsistency。

## 17. 最短复述

Network beam 的算法可以这样复述：

```text
先把每个样本的表达变成 logCPM。
再用训练好的 gene-wise projector 把 logCPM 映射到 Bo2023-like VSD。
然后把 query 和每个 Network 的 centroid 做 Pearson correlation。
分数最高的三个 Network 就是 Top3 beam。
验证时 LOSO/LOMO 都严格在训练折内建立 reference，避免数据泄漏。
这个 beam 是后续 logCPM local reranking 的候选范围，不等于 exact-region 定位结论。
```
