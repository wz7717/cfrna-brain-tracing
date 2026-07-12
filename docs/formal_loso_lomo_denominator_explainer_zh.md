# Formal LOSO/LOMO 分母与 unsupported samples 解释

本文档专门解释 cfRNA-BrainTrace v0.1.7 中最容易混淆的一个问题：

```text
为什么 Bo2023 formal LOSO/LOMO 的 Network 分母是 819，
但 resolution-group / exact-region 分母分别是 814 和 812？
```

如果只记一句话：

```text
每个样本都有 Network truth label，所以 Network 可以评价全部 819 个样本；但 region 层必须要求 held-out 样本的真实 region 在训练 reference 中仍然存在，缺少合法 region reference 的样本不能算对，也不能算错，只能从 group/exact 分母中排除并披露。
```

## 1. LOSO 和 LOMO 是什么

### 1.1 LOSO

LOSO = leave one sample out。

每次只留出 1 个样本做测试，其余样本做训练。

对 819 个 Bo2023 样本来说：

```text
第 1 折:
  test = sample_1
  train = other 818 samples

第 2 折:
  test = sample_2
  train = other 818 samples

...

一共 819 折
```

LOSO 主要测试：

> 如果只拿掉一个样本，模型能不能从其余 reference 中找回它的标签？

### 1.2 LOMO

LOMO = leave one monkey out。

每次留出一整只 monkey 的全部样本，其余 monkeys 做训练。

例如：

```text
第 1 折:
  test = monkey_A 的所有样本
  train = 其他 monkeys 的样本

第 2 折:
  test = monkey_B 的所有样本
  train = 其他 monkeys 的样本
```

LOMO 主要测试：

> 如果训练集中完全没有某只 monkey，模型能不能泛化到这只 monkey？

所以 LOMO 比 LOSO 更严格。

## 2. formal three-tier route 是什么

当前投稿主线是 formal three-tier route：

```text
第 1 层: projected-VSD Network Top3 beam
第 2 层: logCPM-compatible resolution-group reranking
第 3 层: logCPM-compatible exact-region reranking
```

对应脚本：

```text
LOSO:
  scripts/run_bo2023_hybrid_formal_loso.py

LOMO:
  scripts/run_bo2023_projected_vsd_formal_lomo.py
```

这两个脚本都会输出三层结果：

```text
Network
Resolution group
Exact region
```

但这三层的“可评价条件”不完全一样。

## 3. 为什么 Network 可以评价全部 819 个样本

Network 是粗粒度标签。Bo2023 的 819 个样本都有 `SaleemNetworks` truth label。

也就是说，每个 held-out sample 都能回答这个问题：

```text
真实 Network 是否在预测 Top3 里？
```

即使这个样本的 exact region 在训练折里不存在，Network 仍然可以评价，因为 Network reference 仍然存在。

所以：

```text
Formal LOSO Network n = 819
Formal LOMO Network n = 819
```

代码层面也体现了这个逻辑。

在 LOSO 脚本中，Network 结果先写入：

```text
network_rows.append(network_row(...))
```

之后才检查 region 是否可评价。

这意味着：

```text
region 不可评价，不会影响 Network 分母。
```

## 4. 为什么 region 层不能评价全部 819 个样本

region 层需要更细的 truth label，例如具体 Bo2023 region。

要评价一个 held-out sample 的 region prediction，训练集中必须有这个 truth region 的 reference。

否则会出现一个根本问题：

> 模型的候选 reference 里没有这个真实 region，怎么判断它预测错了还是对了？

如果真实 region 不在训练集中，那么算法不可能把它排出来。把这种情况算作错误会惩罚一个“无可评价 reference”的样本；把它算作正确更不合理。

因此唯一严谨做法是：

```text
Network 仍然评价；
region-level group/exact 不评价；
把样本写入 unsupported samples 表；
在论文和 workflow 中披露分母。
```

## 5. LOSO 为什么是 814 个 reference-supported 样本

LOSO 中有 819 个样本。

其中 5 个样本在 leave-one-sample-out 后，训练折中没有对应 truth region reference。

因此：

```text
Network:
  819 个样本全部评价

Resolution group:
  814 个 reference-supported 样本评价

Exact region:
  814 个 reference-supported 样本评价

Unsupported region samples:
  5 个
```

论文 Table S1 写法：

```text
Network n=819; group/exact n=814 reference-supported samples
```

Table S2 当前数字：

```text
Formal LOSO Network:
  n = 819
  Top1 = 58.24%
  Top3 = 92.19%

Formal LOSO Resolution group:
  n = 814
  Top1 = 44.47%
  Top3 = 72.36%

Formal LOSO Exact region:
  n = 814
  Top1 = 22.48%
  Top3 = 45.33%
```

## 6. LOMO 为什么是 812 个 reference-supported 样本

LOMO 中仍然一共有 819 个样本参与 Network 评价。

但由于每次留出的是整只 monkey，训练集中缺失某些 truth region 的情况更容易出现。

有 7 个样本的 truth region 在所有 training monkeys 中不存在。

因此：

```text
Network:
  819 个样本全部评价

Resolution group:
  812 个 reference-supported 样本评价

Exact region:
  812 个 reference-supported 样本评价

Unsupported region samples:
  7 个
```

Table S1 写法：

```text
Network n=819; group/exact n=812 reference-supported samples
```

Table S2 当前数字：

```text
Formal LOMO Network:
  n = 819
  Top1 = 57.75%
  Top3 = 91.21%

Formal LOMO Resolution group:
  n = 812
  Top1 = 41.38%
  Top3 = 69.09%

Formal LOMO Exact region:
  n = 812
  Top1 = 22.17%
  Top3 = 42.36%
```

## 7. unsupported sample 到底是什么

unsupported sample 不是坏样本，也不是预测失败样本。

它的意思是：

```text
这个 held-out 样本的真实 region，在当前训练折的 reference 中不存在。
```

更具体地说：

- LOSO 中：拿掉这个样本后，训练集中没有同 region 的其他样本。
- LOMO 中：拿掉整只 monkey 后，其他 monkeys 中没有这个 region。

因此它是“评价设计限制”，不是“模型错误”。

## 8. 为什么不能把 unsupported sample 算错

假设一个样本真实 region 是 `Region_X`。

在某个 fold 里，训练 reference 中没有任何 `Region_X`。

算法候选列表只能从训练 reference 里的 regions 产生：

```text
候选 regions = [Region_A, Region_B, Region_C, ...]
```

如果 `Region_X` 不在候选空间里，算法不可能预测 `Region_X`。

这时如果强行算错，就等于问：

> 你为什么没有预测一个训练 reference 中不存在的标签？

这是不公平的，也不是严格验证。

## 9. 为什么也不能把 unsupported sample 算对

同样，因为没有合法 reference，算法也没有真正完成 region-level recovery。

所以不能因为它 Network 对了，就把 region 也算对。

正确处理是第三种：

```text
region-level not evaluable
```

也就是：

```text
不进 region 分母；
单独输出 unsupported list；
论文披露数量和原因。
```

## 10. 代码怎样实现这个规则

### 10.1 LOSO 脚本

文件：

```text
scripts/run_bo2023_hybrid_formal_loso.py
```

关键变量：

```text
network_rows
exact_rows
group_rows
unsupported_rows
```

核心顺序：

```text
1. 计算 Network Top3
2. 写入 network_rows
3. 检查 truth_region 是否在 region_training
4. 如果不在:
     写入 unsupported_rows
     continue
5. 如果在:
     写入 exact_rows 和 group_rows
```

这保证：

```text
len(network_rows) = 819
len(group_rows) = 814
len(exact_rows) = 814
len(unsupported_rows) = 5
```

输出文件：

```text
hybrid_formal_loso_region_unsupported_samples.csv
hybrid_formal_loso_summary.json
```

`summary.json` 里会写：

```text
evaluation_denominators:
  network_n
  resolution_group_n
  exact_region_n
  region_unsupported_n
  region_unsupported_reason
```

### 10.2 LOMO 脚本

文件：

```text
scripts/run_bo2023_projected_vsd_formal_lomo.py
```

LOMO 的 unsupported 判断发生在每个 held-out monkey fold 中：

```text
if truth_region not in region_training:
    unsupported_rows.append(...)
```

输出文件：

```text
formal_lomo_region_unsupported_samples.csv
formal_lomo_validation_summary.json
```

`summary.json` 中同样记录：

```text
network_n
resolution_group_n
exact_region_n
region_unsupported_n
region_unsupported_reason
```

LOMO 的 reason 是：

```text
truth region absent from all training monkeys
```

## 11. 旧 92.38% 为什么不能再用

旧的 `92.38%` 来自一个不一致分母：

```text
Formal LOSO Network Top3 只在 814 个 region-evaluable 样本上计算
```

这会造成问题：

```text
Network 是 814 分母
group/exact 也是 814 分母
```

看起来整齐，但不严谨，因为 Network 明明可以评价全部 819 个样本。

当前修正后的规则是：

```text
Network 用 819
group/exact 用 814
```

所以当前投稿数字是：

```text
Formal LOSO Network Top3 = 92.19%
```

而不是：

```text
92.38%
```

`92.38%` 只能作为 legacy denominator inconsistency 披露。

## 12. 为什么分母不同不是“选择性报告”

分母不同是因为任务层级不同：

```text
Network task:
  truth label exists for all 819 samples

Region task:
  truth region must also exist in training reference
```

这不是为了让数字好看，而是为了避免错误定义评价对象。

事实上，分母披露让结果更保守、更透明：

- 没有把 unsupported region 样本伪装成可评价。
- 没有把它们算作成功。
- 也没有把它们算作失败。
- 单独输出文件，让审稿人可以检查。

## 13. 怎么读 summary 文件

### 13.1 LOSO summary

文件：

```text
hybrid_formal_loso_summary.json
```

重点看：

```text
evaluation_denominators.network_n
evaluation_denominators.resolution_group_n
evaluation_denominators.exact_region_n
evaluation_denominators.region_unsupported_n
network.top3_accuracy
resolution_group.group_top3_accuracy
exact_region.top3_accuracy
```

期望：

```text
network_n = 819
resolution_group_n = 814
exact_region_n = 814
region_unsupported_n = 5
```

### 13.2 LOMO summary

文件：

```text
formal_lomo_validation_summary.json
```

重点看：

```text
evaluation_denominators.network_n
evaluation_denominators.resolution_group_n
evaluation_denominators.exact_region_n
evaluation_denominators.region_unsupported_n
```

期望：

```text
network_n = 819
resolution_group_n = 812
exact_region_n = 812
region_unsupported_n = 7
```

## 14. 怎么读 detail 表

### 14.1 Network detail

LOSO：

```text
hybrid_formal_loso_network_detail.csv
```

LOMO：

```text
formal_lomo_network_detail.csv
```

应当有 819 行当前主 route 样本评价。

关键列：

```text
sample_id
label / truth_network
pred_top1
pred_top2
pred_top3
hit1
hit3
true_rank
```

### 14.2 region detail

LOSO：

```text
hybrid_formal_loso_resolution_group_detail.csv
hybrid_formal_loso_exact_region_detail.csv
```

LOMO：

```text
formal_lomo_resolution_group_detail.csv
formal_lomo_exact_region_detail.csv
```

这些表只包含 reference-supported 样本：

```text
LOSO = 814 行
LOMO = 812 行
```

### 14.3 unsupported list

LOSO：

```text
hybrid_formal_loso_region_unsupported_samples.csv
```

LOMO：

```text
formal_lomo_region_unsupported_samples.csv
```

关键列：

```text
sample_id
truth_network
truth_region
reason
network_included_in_evaluation
resolution_group_included_in_evaluation
exact_region_included_in_evaluation
```

期望逻辑：

```text
network_included_in_evaluation = True
resolution_group_included_in_evaluation = False
exact_region_included_in_evaluation = False
```

## 15. LOSO 与 LOMO 为什么不完全同构

LOSO 和 LOMO 都是留出验证，但难度和训练折不同。

LOSO：

```text
只留出一个样本。
训练集中通常还有同一 monkey 的其他样本。
```

LOMO：

```text
留出整只 monkey。
训练集中没有该 monkey 的任何样本。
```

因此 LOMO 更考验 donor-level generalization。

这也是为什么 LOMO 的 region unsupported 样本数是 7，而 LOSO 是 5。

## 16. 审稿复现时的检查清单

### 16.1 检查 Table S1

应看到：

```text
Formal three-tier LOSO:
  Network n=819; group/exact n=814 reference-supported samples

Formal three-tier LOMO:
  Network n=819; group/exact n=812 reference-supported samples
```

### 16.2 检查 Table S2

应看到：

```text
LOSO Network:
  n = 819
  Top1/Top3 = 58.24% / 92.19%

LOSO Resolution group:
  n = 814
  Top1/Top3 = 44.47% / 72.36%

LOSO Exact:
  n = 814
  Top1/Top3 = 22.48% / 45.33%

LOMO Network:
  n = 819
  Top1/Top3 = 57.75% / 91.21%

LOMO Resolution group:
  n = 812
  Top1/Top3 = 41.38% / 69.09%

LOMO Exact:
  n = 812
  Top1/Top3 = 22.17% / 42.36%
```

### 16.3 检查 workflow summary

最终 summary 应明确写：

```text
target
actual
consistent
difference_reason
```

且不应把 `92.38%` 当作当前 route。

### 16.4 检查 unsupported samples 文件

确认：

```text
LOSO unsupported = 5
LOMO unsupported = 7
```

这些样本应被解释为：

```text
region-level not evaluable
```

而不是：

```text
wrong prediction
```

## 17. 最短复述

Formal LOSO/LOMO 的分母规则可以这样复述：

```text
Network 是粗标签，819 个样本都有 truth，所以 Network 全部评价。
Region 是细标签，必须要求 truth region 在训练 fold 中有 reference。
LOSO 有 5 个样本缺 region reference，所以 group/exact 用 814。
LOMO 有 7 个样本缺 region reference，所以 group/exact 用 812。
这些 unsupported 样本不算对、不算错，而是单独披露。
旧 92.38% 是把 Network 也限制到 814 分母的 legacy inconsistency；当前投稿 Network LOSO Top3 是 92.19%。
```
