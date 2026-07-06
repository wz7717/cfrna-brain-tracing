# cfRNA-BrainTrace v0.1.6 逐脚本代码原理讲解

本文档是 `docs/validation_and_tracing_principles_zh.md` 的代码阅读版。前一份文档讲“路线和原理”，本文按脚本解释“代码具体怎样把输入变成论文数字”。

阅读顺序建议：

1. 先看 `core/reference_projection.py`，理解表达矩阵和 projector。
2. 再看 `core/network_tracing.py`，理解 Network beam。
3. 再看 `core/bo2023_region_tracing.py`，理解 region local reranking。
4. 最后看各个 `scripts/run_*` 或 `scripts/analyze_*` 验证脚本。

## 1. 代码里的基本数据形态

### 1.1 expression matrix

大多数脚本使用的表达矩阵形态是：

```text
rows    = genes
columns = samples
values  = counts / logCPM / VSD / projected-VSD
```

也就是一个 `pandas.DataFrame`：

```text
              sample_1  sample_2  sample_3
gene_symbol
GAD1             ...
SLC17A7          ...
MBP              ...
```

### 1.2 metadata

Bo2023 样本 metadata 至少要给出：

```text
sample_id
region_id
network_id
monkey_id   # LOMO 需要
```

在 Excel 里通常对应：

```text
No.
Region
SaleemNetworks
MonkeyID
```

### 1.3 ranking row

验证脚本最后常把每个样本变成一行：

```text
sample_id
label / truth_network / truth_region
pred_top1
pred_top2
pred_top3
hit1
hit3
true_rank
```

其中：

- `hit1 = 1` 表示 Top1 命中。
- `hit3 = 1` 表示真实标签在前三名里。
- `true_rank` 是真实标签在完整排序中的位置。

### 1.4 summary row

汇总时基本就是：

```python
top1_accuracy = detail["hit1"].mean()
top3_accuracy = detail["hit3"].mean()
median_true_rank = detail["true_rank"].median()
```

所以 Top3 的分母就是 detail 表的行数。这也是为什么 Network 和 region 分母必须分开。

## 2. `core/reference_projection.py`

### 2.1 这个模块解决什么问题

这个模块负责把 Bo2023 原始表达输入整理成算法可用的矩阵，并提供 logCPM-to-VSD projector。

核心问题是：

> 用户或外部数据通常只有 raw counts / logCPM，但 Bo2023 reference 的 Network model 是在 VSD-like 空间里效果更好。怎样把 logCPM 映射到 VSD-like 空间？

答案是对每个 gene 拟合一个简单线性映射：

```text
projected_VSD_gene = slope_gene * logCPM_gene + intercept_gene
```

### 2.2 `read_bo2023_gene_matrix`

作用：读取 Bo2023 gene x sample 矩阵。

它处理一个特殊格式：文件第一行是 sample ids，但没有 gene id 列名。因此代码手动加：

```python
names = ["gene_id", *samples]
```

然后读成：

```text
index   = gene_id
columns = sample ids
```

### 2.3 `read_gene_map`

作用：读取 gene id 到 gene symbol 的表。

关键点：

- 读取 `Gene.stable.ID` 和 `Gene.name`。
- 把列名改为 `gene_id` 和 `gene_symbol`。
- 调用 `clean_excel_date_gene_symbol` 修复 Excel 把基因名当日期的问题。

例如有些基因可能被 Excel 误写成类似日期的字符串，代码会把它修回 MARCH/SEPT/DEC 类 symbol。

### 2.4 `map_index_to_symbols`

作用：把 expression matrix 的 index 从 gene id 改成 gene symbol。

如果多个 gene id 映射到同一个 symbol，就按 symbol 聚合取平均：

```python
mapped = matrix.groupby(symbols, sort=True).mean()
```

这一步很重要，因为后面的模型和外部输入通常按 gene symbol 对齐。

### 2.5 `compute_logcpm`

作用：把 raw counts 转成 logCPM。

公式：

```text
library_size = 每个样本所有 gene counts 之和
CPM = counts / library_size * 1,000,000
logCPM = log1p(CPM)
```

代码里还处理了 library size 为 0 的情况，避免除零。

### 2.6 `align_matrices`

作用：让两个矩阵对齐到共同 genes 和共同 samples。

例如 counts 和 VSD 要同时用于训练 projector，必须保证：

```text
same genes
same samples
same order
```

否则每个 gene 的 logCPM 和 VSD 就会错配。

### 2.7 `fit_linear_projector`

作用：为每个 gene 拟合 logCPM 到 VSD 的线性关系。

对每个 gene，代码计算：

- logCPM 均值和方差
- VSD 均值和方差
- slope
- intercept
- residual
- r2
- Spearman correlation
- clip_low / clip_high

如果某个 gene 低表达或 logCPM 方差太低，就不强行拟合 slope，而是 fallback 到 VSD 平均值。

这样做的目的：避免极低信息量 gene 的线性映射产生不稳定外推。

### 2.8 `apply_projector`

作用：把 query logCPM 映射到 projected-VSD。

过程：

1. query 矩阵按 projector genes 重排。
2. 缺失 gene 填 0。
3. 应用每个 gene 的 slope 和 intercept。
4. 用训练集 VSD 分位数裁剪，避免离谱外推。

### 2.9 `save_projector_npz` / `load_projector_npz`

projector 保存为：

```text
data/models/bo2023_reference_projector_linear_full.npz
```

里面包含：

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
metadata
```

## 3. `core/network_tracing.py`

### 3.1 这个模块解决什么问题

生产环境中，用户上传一个表达表，系统要输出 SaleemNetworks 排名。

核心流程是：

```text
user expression
    -> logCPM-like series
    -> projected-VSD
    -> correlate with Network centroids
    -> optional pairwise rescue
    -> Network ranking
```

### 3.2 `load_network_model`

读取：

```text
data/models/bo2023_saleem_network_top200_model.npz
```

模型里有：

```text
genes
networks
reference
fisher_scores
```

这里的 `reference` 可以理解为每个 Network 的 centroid 表达特征。

### 3.3 `_sample_logcpm_series`

作用：把用户上传表格转成一个 gene -> value 的序列。

优先级：

1. 如果有 `read_count`
   - 转成 logCPM。
   - 标记 `input_scale = read_count_logcpm`。

2. 如果有 `log_tpm`
   - 作为 fallback。
   - 标记 `stored_log_tpm_fallback`。

3. 否则用 `tpm_value`
   - 做 `log1p(tpm)`。
   - 标记 `log1p_tpm_fallback`。

投稿主线最推荐 raw counts/logCPM。TPM/logTPM 是兼容旧输入，不应当解释成 Bo2023 当前验证路线完全等价。

### 3.4 `_load_reference_projector`

读取：

```text
data/models/bo2023_reference_projector_linear_full.npz
```

如果文件存在且 `project_to_vsd=True`，Network scoring 使用 projected-VSD。

### 3.5 `trace_network_expression`

这是 Network tracing 主入口。

代码逻辑：

1. 加载 Network model。
2. 把用户 expression 转成 logCPM-like `input_series`。
3. 加载 projector。
4. 如果启用 projector：
   - `apply_projector(projector, query)`
   - 得到 projected-VSD series。
5. 检查 model genes overlap。
6. 如果 overlap 不够，返回 `traceability = insufficient`。
7. 如果 overlap 足够：
   - query vector 按 model genes 重排。
   - 与每个 Network reference 做 Pearson correlation。
   - 分数越高，排名越靠前。
8. 计算 softmax confidence。
9. 可选 pairwise rescue。
10. 返回 `results` 和 `meta`。

### 3.6 `_apply_pairwise_rescue`

pairwise rescue 的作用不是重做完整排序，而是在 Top3 内部修正某些训练中常见混淆。

它只比较：

```text
当前 Top1 anchor
vs
Top2/Top3 challenger
```

如果 pair-specific genes 支持 challenger 且 margin 超过阈值，就交换 Top1。

重要边界：

- rescue 不会引入 Top3 外的新 Network。
- 它只改变 Top1 顺序。
- 因此 Network Top3 beam 的解释仍然稳定。

## 4. `core/bo2023_region_tracing.py`

### 4.1 这个模块解决什么问题

Network tracing 只告诉我们粗 Network。这个模块在 Network Top3 beam 内做 secondary region tracing。

输入：

```text
expression
network_output
db_path
atlas_id
```

输出：

```text
region candidates
resolution groups
scores
traceability metadata
```

### 4.2 `_load_raw_logcpm_reference_matrix`

从 controlled Bo2023 raw counts 加载 reference：

1. 读 counts。
2. 读 sample metadata。
3. 读 gene map。
4. gene id 映射到 symbol。
5. 计算 logCPM。
6. 得到 region/reference matrix。

这就是为什么 region 层是 logCPM-compatible，而不是 projected-VSD direct exact scoring。

### 4.3 `_load_db_reference_matrix`

如果 raw Bo2023 文件不可用，生产 app 可从 SQLite atlas 读取 legacy/reference matrix。

这更多是 app 兼容路径；投稿验证关注 controlled raw logCPM route。

### 4.4 `_candidate_gene_order`

作用：给 candidate regions 选局部有区分度的 genes。

直观理解：

> 如果只在 Network Top3 beam 内比较几个候选脑区，就不需要全基因；更需要那些能区分这些候选脑区的 genes。

### 4.5 `_zscore`

把不同 gene set 的 correlation scores 标准化，使 top50 和 top100 分数组合时尺度可比。

### 4.6 `trace_bo2023_secondary_regions`

主入口逻辑：

1. 从 `network_output["results"]` 取 Network Top3。
2. 如果没有 Network beam，返回空结果和 `traceability = insufficient`。
3. 加载 Bo2023 reference matrix。
4. 筛出 Network Top3 内的 candidate regions。
5. 把 query expression 转成 logCPM-like series。
6. 按 reference genes 对齐 query。
7. 选 local genes。
8. 分别计算 top50 和 top100 correlation。
9. 对两个 score 各自 z-score。
10. 用 fusion weight 合并。
11. 输出 ranked region results。
12. 同时生成 resolution group ordering。

默认 exact route 名：

```text
top3_beam_local_top50_top100_zfusion_w0p25
```

这说明 exact reranking 是：

```text
Network Top3 beam
+ local genes
+ top50/top100 correlation
+ z-score fusion
```

## 5. `scripts/run_bo2023_projected_vsd_loso.py`

### 5.1 这个脚本解决什么问题

验证 projected-VSD 是否适合做 Network beam。

LOSO = leave one sample out。每次留出一个样本，其余样本作为训练集。

### 5.2 入口参数

主要参数：

```text
--counts
--vsd
--sample-info
--sample-sheet
--network-col
--gene-map
--locked-model-genes
--max-folds
--outdir
```

### 5.3 核心函数

`corr_scores(reference, sample)`

- 对 sample 和每个 Network reference 做 Pearson correlation。
- 返回每个 Network 的分数。

`read_labels`

- 从 sample-info Excel 读取 sample_id 和 Network label。

`build_centroids`

- 用训练样本为每个 Network 计算 centroid。

`make_rank_row`

- 把一个 held-out sample 的预测结果写成 detail 表中的一行。

### 5.4 主流程

1. 读 counts、VSD、gene map。
2. gene id 映射到 gene symbol。
3. counts 和 VSD 对齐。
4. counts 转 logCPM。
5. 读取 locked model genes。
6. 对每个 held-out sample：
   - 训练集 = 其他样本。
   - 在训练集内建立 VSD Network centroids。
   - 把 held-out sample 的 logCPM 投影到 VSD-like 空间。
   - 与 centroids 做 correlation。
   - 排序得到 Top1/Top3。
   - 判断真实 Network 是否命中。
7. 汇总所有样本。

### 5.5 输出文件

```text
bo2023_projected_vsd_loso_detail.csv
bo2023_projected_vsd_loso_route_summary.csv
bo2023_projected_vsd_loso_summary.json
```

### 5.6 对应论文数字

```text
n = 819
Top1 = 58.00%
Top3 = 91.58%
median true-rank = 1.0
```

## 6. `scripts/run_bo2023_projected_vsd_lomo.py`

### 6.1 这个脚本解决什么问题

验证 projected-VSD Network beam 是否能跨 monkey 泛化。

LOMO = leave one monkey out。每次留出一只 monkey 的全部样本。

### 6.2 入口参数

主要参数：

```text
--counts
--vsd
--sample-info
--sample-sheet
--monkey-col
--network-col
--gene-map
--locked-model-genes
--outdir
```

### 6.3 主流程

1. 读取 counts、VSD、metadata。
2. 根据 `MonkeyID` 分 fold。
3. 每个 fold：
   - test samples = 当前 monkey 全部样本。
   - train samples = 其他 monkeys。
   - 用 train samples 建立 Network centroids。
   - 对 test samples 做 projected-VSD Network scoring。
   - 记录每个样本的 Top1/Top3。
4. 汇总所有 monkey folds。

### 6.4 输出文件

```text
bo2023_projected_vsd_lomo_detail.csv
bo2023_projected_vsd_lomo_folds.csv
bo2023_projected_vsd_lomo_route_summary.csv
bo2023_projected_vsd_lomo_summary.json
```

### 6.5 对应论文数字

```text
n = 819
Top1 = 53.72%
Top3 = 91.33%
median true-rank = 1.0
```

### 6.6 与 formal LOMO Network 的区别

这个脚本是 Network-only 方法合理性验证。

formal LOMO Network 出自完整三层脚本 `run_bo2023_projected_vsd_formal_lomo.py`。它和 region reranking 绑定在同一 formal route 中，因此作为投稿完整流程的 Network 层结果。

## 7. `scripts/run_bo2023_hybrid_formal_loso.py`

### 7.1 这个脚本解决什么问题

验证完整三层投稿路线在 LOSO 下的表现：

```text
projected-VSD Network beam
-> logCPM resolution-group reranking
-> logCPM exact-region reranking
```

### 7.2 入口参数

主要参数：

```text
--counts
--vsd
--sample-info
--sample-sheet
--region-col
--network-col
--gene-map
--locked-model-genes
--local-top-n-genes
--exact-fusion-weight
--min-resolution-samples
--min-merge-samples
--group-min-pair-errors
--min-confusion-rate
--similarity-threshold
--merge-similarity-threshold
--max-group-size
--max-folds
--outdir
```

### 7.3 核心函数

`read_metadata`

- 读取 sample_id、region_id、network_id。

`network_row`

- 写 Network detail 行。

`exact_row`

- 写 exact-region detail 行。

`summarize_network`

- 统计 Network Top1/Top3。

`summarize_exact`

- 统计 exact-region Top1/Top3。

`summarize_group`

- 统计 resolution-group Top1/Top3，同时保留 exact hit 信息。

### 7.4 主循环的关键逻辑

对每个 held-out sample：

1. 建立训练集。
2. 拿到 truth Network 和 truth region。
3. 用训练集 VSD centroids 做 Network scoring。
4. 立刻把 Network 结果写入 `network_rows`。
5. 检查 truth region 是否还在训练集。
6. 如果 truth region 不在训练集：
   - 写入 `unsupported_rows`。
   - Network 保留。
   - group/exact 跳过。
7. 如果 truth region 在训练集：
   - 在 Network Top3 beam 中找 candidate regions。
   - 选 local discriminative genes。
   - 计算 exact-region top50/top100 z-score fusion。
   - 构建 resolution groups。
   - 计算 group ranking。
   - 写入 exact 和 group detail。

### 7.5 为什么 Network 是 819，而 group/exact 是 814

因为第 4 步在 region support 检查前已经记录 Network。

但有 5 个 held-out samples 的 truth region 在训练折里不存在。它们没有合法 region reference，因此：

```text
Network included = yes
region included = no
```

### 7.6 输出文件

```text
hybrid_formal_loso_network_detail.csv
hybrid_formal_loso_exact_region_detail.csv
hybrid_formal_loso_resolution_group_detail.csv
hybrid_formal_loso_fold_summary.csv
hybrid_formal_loso_region_unsupported_samples.csv
hybrid_formal_loso_network_route_metrics.csv
hybrid_formal_loso_exact_region_route_metrics.csv
hybrid_formal_loso_resolution_group_route_metrics.csv
hybrid_formal_loso_summary.json
```

### 7.7 对应论文数字

```text
Network:
  n = 819
  Top1 = 58.24%
  Top3 = 92.19%

Resolution group:
  n = 814
  Top1 = 44.47%
  Top3 = 72.36%

Exact region:
  n = 814
  Top1 = 22.48%
  Top3 = 45.33%
```

## 8. `scripts/run_bo2023_projected_vsd_formal_lomo.py`

### 8.1 这个脚本解决什么问题

验证完整三层投稿路线在 LOMO 下的表现。

相比 LOSO，它更严格，因为整只 monkey 被留出。

### 8.2 route families

脚本里会同时评估多个 route family：

```text
native_vsd
projected_vsd
logcpm_baseline
hybrid_projected_network_logcpm_exact
```

当前投稿主线关注：

```text
hybrid_projected_network_logcpm_exact
```

### 8.3 主流程

对每个 held-out monkey：

1. `test_indices` = 当前 monkey 的所有样本。
2. `train_indices` = 其他 monkeys 的样本。
3. 建立训练集 region reference。
4. 对 test samples 做 fold-local projection。
5. 对每个 route family：
   - 建立 Network reference。
   - 选 Network discriminative genes。
   - 建立 pairwise rescue models。
   - 对每个 test sample 做 Network scoring。
   - 保存 Network detail。
   - 如果 truth region 不在 training monkeys 中，跳过 region。
   - 如果可评价，做 exact 和 group local reranking。

### 8.4 pairwise rescue 在这里的作用

formal LOMO 的 Network 层包含 pairwise Top1 rescue。它只在 Top3 内部调整排序，不把 Top3 外的新 Network 加进来。

这就是为什么 formal LOMO Network 数字和 independent projected-VSD Network LOMO 数字接近但不完全一样。

### 8.5 输出文件

```text
formal_lomo_fold_summary.csv
formal_lomo_region_unsupported_samples.csv
formal_lomo_network_detail.csv
formal_lomo_exact_region_detail.csv
formal_lomo_resolution_group_detail.csv
formal_lomo_network_route_metrics.csv
formal_lomo_exact_region_route_metrics.csv
formal_lomo_resolution_group_route_metrics.csv
formal_lomo_network_per_monkey_metrics.csv
formal_lomo_exact_region_per_monkey_metrics.csv
formal_lomo_resolution_group_per_monkey_metrics.csv
formal_lomo_validation_summary.json
```

### 8.6 对应论文数字

当前投稿主线 route family：

```text
hybrid_projected_network_logcpm_exact
```

论文数字：

```text
Network:
  n = 819
  Top1 = 57.75%
  Top3 = 91.21%

Resolution group:
  n = 812
  Top1 = 41.38%
  Top3 = 69.09%

Exact region:
  n = 812
  Top1 = 22.17%
  Top3 = 42.36%
```

## 9. `scripts/run_ahba_projected_vsd_formal_three_tier_external.py`

### 9.1 这个脚本解决什么问题

用 AHBA human brain RNA-seq 做 mapped-label external validation。

它回答：

> 当人类 AHBA 解剖标签被映射到当前 macaque-derived hierarchy 后，hybrid route 的候选排序是否与映射标签一致？

### 9.2 入口参数

主要参数包括：

```text
--ahba-expression
--ahba-samples
--bo-counts
--bo-vsd
--sample-info
--region-col
--network-col
--gene-map
--projector
--global-top-n-genes
--network-gene-pool-size
--local-top-n-genes
--exact-fusion-weight
--outdir
```

### 9.3 主流程

1. 读取 AHBA expression 和 sample labels。
2. 读取 Bo2023 counts/VSD/reference metadata。
3. 建立 Bo2023 Network 和 region references。
4. 将 AHBA genes 与 Bo2023 model genes 对齐。
5. 对 AHBA query 做 projected-VSD Network scoring。
6. 对 mapped-label 支持的样本做 group/exact reranking。
7. 输出 sample detail、metrics 和特殊标签审计。

### 9.4 输出文件

```text
ahba_formal_three_tier_sample_detail.csv
ahba_formal_three_tier_metrics.csv
ahba_formal_three_tier_special_labels.csv
ahba_formal_three_tier_resolution_audit.csv
ahba_formal_three_tier_summary.json
```

### 9.5 对应论文数字

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
```

解释边界：这是 mapped-label transfer，不是 human-to-macaque exact anatomical identity。

## 10. `scripts/run_tcga_labeled_hybrid_formal_external.py`

### 10.1 这个脚本解决什么问题

用 TCGA/BraTS glioma tissue RNA-seq 和 MRI-derived labels 做 coarse consistency。

它不证明 macaque exact-region localization，只检查输出是否在粗解剖层面与 MRI label 有一致性。

### 10.2 入口参数

主要参数：

```text
--tcga-counts
--manifest
--labels
--bo-counts
--bo-vsd
--sample-info
--region-col
--network-col
--gene-map
--projector
--region-meta
--global-top-n-genes
--network-gene-pool-size
--local-top-n-genes
--exact-fusion-weight
--outdir
```

### 10.3 核心函数

`split_candidates`

- 把 truth label 中可能的多个候选拆开。

`hit`

- 判断预测列表与 truth candidates 是否有交集。

`load_tcga_counts`

- 读取 TCGA count matrix。

`candidate_regions`

- 根据 Network Top3 beam 选 Bo2023 candidate regions。

`summarize`

- 汇总 Network/lobe/broad anatomy Top1/Top3。

### 10.4 主流程

1. 读取 TCGA counts、manifest、MRI labels。
2. 读取 Bo2023 reference。
3. 计算 TCGA query 的 logCPM。
4. 用 projector 生成 Network Top3 beam。
5. 在 Network beam 内做 local region/group ranking。
6. 将预测结果映射到 human-comparable coarse labels：
   - Network
   - lobe
   - broad anatomy
7. 用 `hit` 判断 Top1/Top3 是否与 truth candidates 相交。
8. 输出 detail、metrics、summary 和 plot。

### 10.5 输出文件

```text
tcga_labeled_hybrid_formal_sample_detail.csv
tcga_labeled_hybrid_formal_metrics.csv
tcga_labeled_hybrid_formal_summary.json
tcga_labeled_hybrid_formal_metrics.png
```

### 10.6 对应论文数字

当前 Table S4 route：

```text
n = 65
Network Top1 = 15.38%
Network Top3 = 40.00%
Broad anatomy Top1 = 13.85%
Broad anatomy Top3 = 64.62%
```

旧 side route：

```text
Network Top3 = 36.92%
Broad Top3 = 80.00%
```

旧 side route 不能作为当前投稿 route。

## 11. `scripts/analyze_gse189919_csf_tracing_validation.py`

### 11.1 这个脚本解决什么问题

GSE189919 没有 anatomical truth，所以它不能计算 localization accuracy。

这个脚本只回答：

> 外部 CSF expression matrix 是否有足够 gene overlap，使当前 projection/tracing pipeline 技术上可以运行和审计？

### 11.2 输入

默认读取：

```text
GSE189919_count.csv.gz
GSE189919_tpm_count.csv.gz
soft metadata
```

入口参数：

```text
--data-dir
--outdir
```

### 11.3 核心函数

`read_soft_metadata`

- 读取样本 metadata。

`read_expression`

- 读取 counts 和 TPM 矩阵。
- 对 sample id 做标准化。

`production_predictions`

- 用 production Network tracing 生成 Network Top1/Top3。

`adapted_production_predictions`

- 对不同输入尺度或策略做适配预测。

`sample_detection`

- 计算 sample-level input QC。

`prediction_distribution`

- 统计各 route 输出的 Network/broad category 分布。

`route_agreement`

- 比较不同 route 的 Top3 Jaccard overlap。

`audit_algorithms`

- 检查旧 baseline/production route 的一些技术假设。

`methodological_review`

- 输出方法学解释边界。

### 11.4 主流程

1. 读 counts、TPM、metadata。
2. 对 sample 和 gene 进行清理对齐。
3. 计算 projector/model gene overlap。
4. 运行若干 prediction routes。
5. 输出预测分布和 route agreement。
6. 运行 optional algorithm audit。
7. 输出 methodological review。

### 11.5 输出文件

```text
gse189919_csf_network_predictions.csv
sample_input_qc.csv
sample_detection_summary.csv
network_broad_top1_occurrence.csv
network_broad_top3_occurrence.csv
network_broad_prediction_collapse.csv
route_agreement.csv
algorithm_audit.csv
methodological_review.csv
```

### 11.6 对应论文数字

```text
Projector gene overlap = 15622 / 21668
Coverage = 72.10%
Conclusion = projection feasibility only
```

### 11.7 algorithm audit 失败的正确解释

如果 optional audit 报：

```text
Algorithm audit failed:
baseline_production_tpm ...
baseline_production_count_cpm ...
```

这不是当前投稿准确率路线失败。

它表示旧 baseline/production comparison 的假设不应当作为当前 localization accuracy 证据。workflow 应把它记为：

```text
EXPECTED_LEGACY_ALGORITHM_LIMITATION
```

## 12. `scripts/stage_validation_data_under_tests.py`

### 12.1 这个脚本解决什么问题

公开 GitHub 不能包含 Bo2023 controlled expression matrix。Victor 服务器上如果有授权数据，需要把它们集中到项目下：

```text
tests/controlled_data/
```

这个脚本就是做 controlled data staging。

### 12.2 安全边界

它只允许从 `/storage/wangzhen` 范围内查找和复制文件，避免 workflow 为了复现访问其他用户目录。

### 12.3 典型目标结构

```text
tests/controlled_data/bo2023/
tests/controlled_data/ahba/
tests/controlled_data/tcga_brats/
tests/controlled_data/gse189919/
```

### 12.4 与 pytest 的关系

如果 controlled Bo2023 文件存在：

```text
17 passed
```

如果公开 clone 缺少 controlled Bo2023 文件：

```text
16 passed, 1 skipped
```

这两种结果都合理。

## 13. HTML workflow 与这些脚本的关系

HTML workflow 文件：

```text
output/cfrna_braintrace_v016_jupyter_validation_workflow.html
```

它不是新的算法实现，而是把上面的脚本按 Jupyter Lab cell 顺序串起来：

```text
环境检查
pytest tests
projected-VSD Network LOSO
projected-VSD Network LOMO
formal three-tier LOSO
formal three-tier LOMO
AHBA mapped-label external validation
TCGA/BraTS coarse consistency
GSE189919 projection feasibility
legacy/inconsistency check
final summary
```

最终输出：

```text
output/v016_jupyter_validation/jupyter_validation_summary.csv
output/v016_jupyter_validation/jupyter_validation_summary.md
```

## 14. 从代码到论文数字的对应关系

| 论文结果 | 脚本 | 关键输出 |
|---|---|---|
| Bo2023 projected-VSD Network LOSO 58.00%/91.58% | `run_bo2023_projected_vsd_loso.py` | `bo2023_projected_vsd_loso_route_summary.csv` |
| Bo2023 projected-VSD Network LOMO 53.72%/91.33% | `run_bo2023_projected_vsd_lomo.py` | `bo2023_projected_vsd_lomo_route_summary.csv` |
| Formal LOSO Network 58.24%/92.19% | `run_bo2023_hybrid_formal_loso.py` | `hybrid_formal_loso_network_route_metrics.csv` |
| Formal LOSO group 44.47%/72.36% | `run_bo2023_hybrid_formal_loso.py` | `hybrid_formal_loso_resolution_group_route_metrics.csv` |
| Formal LOSO exact 22.48%/45.33% | `run_bo2023_hybrid_formal_loso.py` | `hybrid_formal_loso_exact_region_route_metrics.csv` |
| Formal LOMO Network 57.75%/91.21% | `run_bo2023_projected_vsd_formal_lomo.py` | `formal_lomo_network_route_metrics.csv` |
| Formal LOMO group 41.38%/69.09% | `run_bo2023_projected_vsd_formal_lomo.py` | `formal_lomo_resolution_group_route_metrics.csv` |
| Formal LOMO exact 22.17%/42.36% | `run_bo2023_projected_vsd_formal_lomo.py` | `formal_lomo_exact_region_route_metrics.csv` |
| AHBA Network 74.68%/94.42% | `run_ahba_projected_vsd_formal_three_tier_external.py` | `ahba_formal_three_tier_metrics.csv` |
| AHBA group 36.26%/67.03% | `run_ahba_projected_vsd_formal_three_tier_external.py` | `ahba_formal_three_tier_metrics.csv` |
| AHBA exact 24.18%/42.86% | `run_ahba_projected_vsd_formal_three_tier_external.py` | `ahba_formal_three_tier_metrics.csv` |
| TCGA/BraTS Network 15.38%/40.00% | `run_tcga_labeled_hybrid_formal_external.py` | `tcga_labeled_hybrid_formal_metrics.csv` |
| TCGA/BraTS broad 13.85%/64.62% | `run_tcga_labeled_hybrid_formal_external.py` | `tcga_labeled_hybrid_formal_metrics.csv` |
| GSE189919 overlap 15622/21668, 72.10% | `analyze_gse189919_csf_tracing_validation.py` | workflow display table / GSE summary outputs |

## 15. 读代码时最重要的三个判断

### 15.1 这个结果是不是当前 submission route

看 route family 或脚本文档是否指向：

```text
hybrid_projected_network_logcpm_exact
```

如果是 old baseline、TPM fallback、旧 MRI truth side route，就不能当当前投稿主结果。

### 15.2 这个结果的分母是什么

必须问：

```text
Network denominator?
Region denominator?
Unsupported region samples?
```

Bo2023 当前主线：

```text
LOSO Network = 819
LOSO group/exact = 814
LOMO Network = 819
LOMO group/exact = 812
```

### 15.3 这个数据集有没有 truth

如果没有 truth，就不能算 accuracy。

```text
Bo2023: has macaque internal labels
AHBA: has mapped labels
TCGA/BraTS: has coarse MRI-derived labels
GSE189919: no anatomical truth
```

因此 GSE189919 只能是 projection feasibility。

## 16. 最短理解版本

如果只记一条线：

```text
raw counts
-> logCPM
-> projected-VSD only for Network Top3
-> logCPM local reranking for group/exact
-> Network metrics use all valid Network labels
-> region metrics only use reference-supported samples
-> external datasets only claim what their truth labels support
```

这就是 cfRNA-BrainTrace v0.1.6 当前投稿验证流程和溯源算法的核心。
