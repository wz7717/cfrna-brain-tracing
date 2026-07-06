# cfRNA-BrainTrace v0.1.6 验证流程与溯源算法原理说明

本文档面向没有读过项目代码的审稿人或合作者，解释 cfRNA-BrainTrace v0.1.6 的溯源算法、验证路线、脚本输入输出、论文数字来源和常见误解边界。

核心结论先写在前面：

- 当前投稿主线不是旧 baseline。当前主线是 projected-VSD 生成 Network Top3 beam，再用 logCPM-compatible 表达做 resolution-group 和 exact-region local reranking。
- projected-VSD 只用于粗粒度 Network 候选束。下游 region 层不把 projected-VSD 当成直接 exact-region accuracy endpoint。
- Network 指标使用全部 Network-evaluable 样本。Bo2023 内部验证中为 819 个样本。
- region 层指标只评价 reference-supported 样本。LOSO 中为 814 个样本，LOMO 中为 812 个样本。
- GSE189919 没有 anatomical truth，只能验证 projection feasibility，不能作为 localization accuracy。
- TCGA/BraTS 是 coarse consistency，不是 macaque exact-region localization。

## 1. 概念地图

### 1.1 三个解剖层级

项目输出不是单一脑区标签，而是层级式候选排序：

1. Network
   - 粗粒度 SaleemNetworks 类别。
   - 投稿主 endpoint 是 Network Top3。
   - 作用是给后续 local reranking 画出候选范围。

2. Resolution group
   - 在 Network Top3 beam 内，把训练集中容易混淆或分辨率不足的脑区合并成较稳健的局部组。
   - 是更可辩护的 region-level endpoint。

3. Exact region
   - 在 Network Top3 beam 内对具体 Bo2023 region 排名。
   - 当前论文中保留为 exploratory candidate ranking，不作为强定位 endpoint。

### 1.2 两种表达尺度

项目同时用到两种表达表示：

1. logCPM
   - 从 raw counts 计算：
     `CPM = counts / library_size * 1,000,000`
     `logCPM = log1p(CPM)`
   - 适合用户上传 raw counts 或 logCPM 后的可复现处理。
   - 下游 resolution-group 和 exact-region reranking 使用 logCPM-compatible 表达。

2. projected-VSD
   - Bo2023 reference 有 native VSD 表达，但普通外部用户通常没有 VSD。
   - 项目训练了一个 gene-wise linear projector，把 logCPM 近似映射到 Bo2023-like VSD 空间。
   - 代码位置：
     - `core/reference_projection.py`
     - `data/models/bo2023_reference_projector_linear_full.npz`
   - projected-VSD 的用途是生成 Network Top3 beam。

### 1.3 为什么先 Network，再 region

直接在所有 exact regions 上打分会遇到两个问题：

- 很多 region 标签非常细，样本数少，跨 monkey 泛化更难。
- 外部数据常常只有粗标签，不能支持 exact-region accuracy。

因此当前主路线采用：

```text
raw counts / logCPM query
        |
        v
logCPM-to-VSD projector
        |
        v
projected-VSD Network Top3 beam
        |
        v
logCPM-compatible local reranking within that beam
        |
        +--> resolution-group ranking
        +--> exact-region exploratory ranking
```

这就是论文里 "projected-VSD Network beam + logCPM local reranking" 的含义。

## 2. 核心代码模块

### 2.1 `core/reference_projection.py`

作用：处理 Bo2023 expression matrix、gene mapping、logCPM 计算和 logCPM-to-VSD projector。

关键函数：

- `read_bo2023_gene_matrix`
  - 读取 Bo2023 gene x sample 矩阵。

- `read_gene_map`
  - 读取 Ensembl gene id 到 gene symbol 的映射。
  - 修正 Excel 日期误转换造成的基因名问题，例如 MARCH/SEPT/DEC 类基因。

- `map_index_to_symbols`
  - 把矩阵 index 从 gene id 映射到 gene symbol。
  - 重复 symbol 取均值。

- `compute_logcpm`
  - 从 raw counts 得到 logCPM。

- `fit_linear_projector`
  - 对每个 gene 拟合：
    `VSD_gene = slope_gene * logCPM_gene + intercept_gene`
  - 对低表达或方差不足的 gene 使用 fallback。

- `apply_projector`
  - 把 query logCPM 映射到 projected-VSD。

- `load_projector_npz`
  - 加载 `bo2023_reference_projector_linear_full.npz`。

### 2.2 `core/network_tracing.py`

作用：生产环境里的 Network 层溯源。

主要流程：

1. `_sample_logcpm_series`
   - 如果用户上传 `read_count`，内部转为 logCPM。
   - 如果用户上传 `log_tpm` 或 `tpm_value`，作为 legacy fallback，不等同于当前 Bo2023 logCPM route。

2. `trace_network_expression`
   - 读取 Network model：
     `data/models/bo2023_saleem_network_top200_model.npz`
   - 默认加载 projector：
     `data/models/bo2023_reference_projector_linear_full.npz`
   - 把 query logCPM 投影为 projected-VSD。
   - 与每个 Network centroid 做 Pearson correlation。
   - 输出 Network 排名、confidence、overlap、projection metadata。
   - 可使用 pairwise rescue，但 rescue 约束在原始 Top3 候选内部。

重点：这里的 projected-VSD 是 Network 层用的。

### 2.3 `core/bo2023_region_tracing.py`

作用：生产环境里的 Bo2023 secondary region tracing。

主要流程：

1. 接收 Network Top3 beam。
2. 从 Bo2023 raw counts/logCPM reference 中筛选属于这 3 个 Network 的候选 regions。
3. 使用 local discriminative genes。
4. 对候选 regions 做 top50/top100 gene correlation。
5. 使用 z-score fusion 生成 exact-region ranking。

重点：region 层使用 logCPM-compatible reference，不把 projected-VSD 直接当 exact-region 主 endpoint。

### 2.4 `core/region_resolution.py`

作用：支持 resolution-group 的分辨率标注。

核心思想：

- 如果一些 region 在训练集中相似度过高、样本数太少或混淆率高，就不强行把它们拆成 exact region。
- 这些 region 会被合并或标成 low-resolution group。
- 这让 group-level 结果比 exact-region 更稳健。

## 3. 当前投稿验证路线总览

| 路线 | 脚本 | 目的 | 投稿解释 |
|---|---|---|---|
| Bo2023 projected-VSD Network LOSO | `scripts/run_bo2023_projected_vsd_loso.py` | 验证 projected-VSD 能否做 Network beam | 方法合理性验证 |
| Bo2023 projected-VSD Network LOMO | `scripts/run_bo2023_projected_vsd_lomo.py` | 验证 projected-VSD Network beam 跨 monkey 泛化 | 方法合理性验证 |
| Bo2023 formal three-tier LOSO | `scripts/run_bo2023_hybrid_formal_loso.py` | 完整投稿路线的 sample-level 内部验证 | 主内部验证 |
| Bo2023 formal three-tier LOMO | `scripts/run_bo2023_projected_vsd_formal_lomo.py` | 完整投稿路线的 donor-level 内部验证 | 主内部验证 |
| AHBA external | `scripts/run_ahba_projected_vsd_formal_three_tier_external.py` | human mapped-label external validation | mapped-label transfer |
| TCGA/BraTS external | `scripts/run_tcga_labeled_hybrid_formal_external.py` | glioma tissue + MRI label coarse consistency | coarse consistency only |
| GSE189919 | `scripts/analyze_gse189919_csf_tracing_validation.py` | 外部 CSF matrix 投影可行性 | no accuracy claim |
| Jupyter HTML workflow | `output/cfrna_braintrace_v016_jupyter_validation_workflow.html` | 把以上路线串成审稿复现工作流 | public/controlled 双模式 |

## 4. Bo2023 projected-VSD Network LOSO

### 4.1 这条路线问什么问题

问题是：

> 如果每次留出 1 个 Bo2023 样本，用剩余样本训练 Network centroid，那么 projected-VSD query 能否把这个样本的真实 Network 排进 Top3？

它只评价 Network，不评价 resolution group 或 exact region。

### 4.2 主要脚本

`scripts/run_bo2023_projected_vsd_loso.py`

### 4.3 输入

通常需要：

- Bo2023 raw counts
- Bo2023 native VSD matrix
- Bo2023 sample metadata
- cleaned gene map

### 4.4 算法步骤

对每个 held-out sample：

1. 从 raw counts 计算 logCPM。
2. 用训练折拟合或应用 logCPM-to-VSD projection。
3. 用训练折样本为每个 Network 计算 centroid。
4. 对 held-out sample 与每个 Network centroid 做 Pearson correlation。
5. correlation 从高到低排序。
6. 判断真实 Network 是否在 Top1 或 Top3。

### 4.5 输出

典型输出：

- `bo2023_projected_vsd_loso_detail.csv`
  - 每个样本的真实 Network、预测 Top1/Top3、是否命中。

- `bo2023_projected_vsd_loso_route_summary.csv`
  - 汇总 Top1、Top3、median true rank。

- `bo2023_projected_vsd_loso_summary.json`
  - 参数和汇总审计信息。

### 4.6 论文数字

Table S2 当前数字：

```text
n = 819
Top1 = 58.00%
Top3 = 91.58%
median true-rank = 1.0
```

### 4.7 如何解释

这条路线证明 projected-VSD 对粗粒度 Network beam 是有效的。它不是完整三层投稿路线，也不支持 exact-region localization claim。

## 5. Bo2023 projected-VSD Network LOMO

### 5.1 这条路线问什么问题

问题是：

> 如果一次留出一个 monkey 的全部样本，projected-VSD Network beam 是否还能跨 donor 泛化？

它仍然只评价 Network。

### 5.2 主要脚本

`scripts/run_bo2023_projected_vsd_lomo.py`

### 5.3 算法步骤

对每个 held-out monkey：

1. 训练集为其他 monkeys。
2. 测试集为该 monkey 所有样本。
3. 用训练 monkeys 建立 Network centroid。
4. 把测试样本的 logCPM 投影到 projected-VSD。
5. 做 Network correlation ranking。
6. 统计所有 819 个样本的 Top1/Top3。

### 5.4 输出

- `bo2023_projected_vsd_lomo_detail.csv`
- `bo2023_projected_vsd_lomo_folds.csv`
- `bo2023_projected_vsd_lomo_route_summary.csv`
- `bo2023_projected_vsd_lomo_summary.json`

### 5.5 论文数字

```text
n = 819
Top1 = 53.72%
Top3 = 91.33%
median true-rank = 1.0
```

### 5.6 为什么它要和 formal LOMO Network 分开说

这条路线是 projected-VSD Network beam 的独立合理性验证。

formal LOMO Network 是完整三层投稿流程里的 Network 层结果。两者都评价 Network，但上下文不同：

- projected-VSD Network LOMO：只问 projected-VSD Network beam 本身是否成立。
- formal LOMO Network：问完整 three-tier route 在 LOMO 下的 Network 层表现。

因此它们不能混成一个数字。

## 6. Bo2023 formal three-tier LOSO

### 6.1 这条路线问什么问题

问题是：

> 在完整投稿路线中，留出单个样本时，Network、resolution group、exact region 三层分别表现如何？

这是主内部验证之一。

### 6.2 主要脚本

`scripts/run_bo2023_hybrid_formal_loso.py`

### 6.3 关键设计

这条路线有两个分母：

1. Network denominator
   - 全部 819 个样本。
   - 因为每个样本都有 Network truth label。

2. Region denominator
   - 814 个 reference-supported 样本。
   - 有 5 个样本在 LOSO 训练折中没有对应 truth region reference，因此不能评价 group/exact。
   - 它们不算 region 成功，也不算 region 失败。
   - 它们仍然参与 Network 评价。

### 6.4 算法步骤

对每个 held-out sample：

1. 训练集 = 其他 818 个样本。
2. 建立 projected-VSD Network centroids。
3. held-out sample 的 logCPM 投影到 projected-VSD。
4. 与 Network centroids 做 correlation，得到 Network Top3 beam。
5. 立刻记录 Network Top1/Top3 命中。
6. 检查 truth region 是否存在于训练折。
7. 如果不存在：
   - 写入 `hybrid_formal_loso_region_unsupported_samples.csv`
   - 跳过 resolution-group 和 exact-region。
8. 如果存在：
   - 在 Network Top3 beam 对应的 training regions 内构造 candidate regions。
   - 使用 logCPM expression 和 local discriminative genes。
   - exact-region：top50/top100 correlation z-score fusion。
   - resolution-group：根据局部分辨率规则合并或标注 candidates，再计算 group hit。

### 6.5 输出

- `hybrid_formal_loso_network_detail.csv`
  - 819 行 Network 评价。

- `hybrid_formal_loso_network_route_metrics.csv`
  - Network Top1/Top3 汇总。

- `hybrid_formal_loso_resolution_group_detail.csv`
  - 814 行 group 评价。

- `hybrid_formal_loso_resolution_group_route_metrics.csv`
  - group Top1/Top3 汇总。

- `hybrid_formal_loso_exact_region_detail.csv`
  - 814 行 exact 评价。

- `hybrid_formal_loso_exact_region_route_metrics.csv`
  - exact Top1/Top3 汇总。

- `hybrid_formal_loso_region_unsupported_samples.csv`
  - 5 个不能 region-level 评价的样本。

- `hybrid_formal_loso_summary.json`
  - denominator policy、implementation note 和汇总结果。

### 6.6 论文数字

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

### 6.7 旧 92.38% 是什么

旧的 `92.38%` 是把 Network Top3 也条件化到 814 个 region-evaluable 样本上得到的数字。当前投稿不使用它。

正确做法是：

- Network 用 819。
- group/exact 用 814。

所以 `92.38%` 只能作为 legacy denominator inconsistency 披露。

## 7. Bo2023 formal three-tier LOMO

### 7.1 这条路线问什么问题

问题是：

> 在完整投稿路线中，留出整只 monkey 时，三层结果能否跨 donor 泛化？

这是比 LOSO 更严格的内部验证。

### 7.2 主要脚本

`scripts/run_bo2023_projected_vsd_formal_lomo.py`

### 7.3 与 LOSO 的不同

LOSO 留出一个样本，训练集中通常仍有同一 monkey 的其他样本。

LOMO 留出一个 monkey 的所有样本，训练集中没有该 donor 的任何样本。因此 LOMO 更接近跨个体泛化。

### 7.4 关键设计

这条路线也有两个分母：

1. Network denominator
   - 819 个样本。

2. Region denominator
   - 812 个 reference-supported 样本。
   - 7 个样本在所有 training monkeys 中没有对应 truth region reference。

### 7.5 算法步骤

对每个 held-out monkey：

1. 训练集 = 其他 monkeys。
2. 测试集 = held-out monkey 的全部样本。
3. 对测试样本做 fold-local projected-VSD Network scoring。
4. Network 层使用 Top3 pairwise rescue，但 rescue 只在原始 Top3 set 内调整 Top1。
5. 对每个测试样本记录 Network 命中。
6. 如果 truth region 不在 training monkeys 中，跳过 group/exact。
7. 如果 truth region 支持评价：
   - 在 Network Top3 beam 内选 candidate regions。
   - 使用 logCPM-compatible exact values。
   - 计算 exact-region 和 resolution-group ranking。

### 7.6 输出

- `formal_lomo_network_detail.csv`
- `formal_lomo_network_route_metrics.csv`
- `formal_lomo_exact_region_detail.csv`
- `formal_lomo_exact_region_route_metrics.csv`
- `formal_lomo_resolution_group_detail.csv`
- `formal_lomo_resolution_group_route_metrics.csv`
- `formal_lomo_network_per_monkey_metrics.csv`
- `formal_lomo_exact_region_per_monkey_metrics.csv`
- `formal_lomo_resolution_group_per_monkey_metrics.csv`
- `formal_lomo_region_unsupported_samples.csv`
- `formal_lomo_validation_summary.json`

### 7.7 论文数字

```text
Formal hybrid Network:
  n = 819
  Top1 = 57.75%
  Top3 = 91.21%

Formal hybrid Resolution group:
  n = 812
  Top1 = 41.38%
  Top3 = 69.09%

Formal hybrid Exact region:
  n = 812
  Top1 = 22.17%
  Top3 = 42.36%
```

### 7.8 为什么 LOMO Network 和 projected-VSD Network LOMO 数字不同

两者都评价 Network，但路线不同：

- `run_bo2023_projected_vsd_lomo.py` 是独立 Network-only projected-VSD validation。
- `run_bo2023_projected_vsd_formal_lomo.py` 是完整 formal three-tier LOMO，同时保留多个 route families，并使用 formal route 的 Network scoring/rescue 设置。

因此：

```text
projected-VSD Network LOMO Top3 = 91.33%
formal LOMO Network Top3 = 91.21%
```

两者接近，但不能互相替代。

## 8. AHBA mapped-label external validation

### 8.1 这条路线问什么问题

问题是：

> 人类 AHBA RNA-seq 样本的解剖标签映射到 macaque-derived hierarchy 后，当前 hybrid route 能否给出合理候选？

它是外部 mapped-label 验证，不是直接同物种 exact anatomical truth。

### 8.2 主要脚本

`scripts/run_ahba_projected_vsd_formal_three_tier_external.py`

相关辅助脚本：

- `scripts/run_ahba_human_rnaseq_external_validation.py`
- `scripts/run_ahba_projected_vsd_external_validation.py`
- `scripts/export_ahba_external_validation_figures.py`

### 8.3 算法步骤

1. 读取 AHBA human RNA-seq expression。
2. 将 human sample anatomical label 映射到当前 hierarchy 可评价标签。
3. 对 query 表达做与当前 route 兼容的 Network beam generation。
4. 在支持 exact mapping 的样本上做 group/exact reranking。
5. 只在 label mapping 支持的样本上计算 accuracy。
6. 同时报告 baseline route 与 projected-VSD-only exact route 作为对照。

### 8.4 论文数字

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

logCPM baseline Exact Top3:
  30.77%

Projected-VSD-only Exact Top3:
  29.67%
```

### 8.5 如何解释

AHBA 支持的是 mapped-label transfer。由于 human label 被映射到 macaque-derived hierarchy，它证明的是外部标签体系下的候选一致性，不是直接的 macaque exact-region truth。

## 9. TCGA/BraTS coarse consistency

### 9.1 这条路线问什么问题

问题是：

> 对 glioma tissue RNA-seq 样本，在 MRI-derived human brain labels 的粗标签下，模型输出是否有粗解剖一致性？

它不是正常脑组织，也不是 Bo2023 macaque exact-region validation。

### 9.2 主要脚本

`scripts/run_tcga_labeled_hybrid_formal_external.py`

相关脚本：

- `scripts/build_corrected_brats_mri_truth.py`
- `scripts/evaluate_brats_tcga_lgg_65_mri_truth.py`
- `scripts/analyze_tcga_gbm_lgg_network_domain_shift.py`

### 9.3 算法步骤

1. 读取 TCGA/BraTS expression。
2. 读取或构建 MRI-derived truth labels。
3. 对样本做 hybrid route prediction：
   - projected-VSD Network Top3 beam
   - logCPM local reranking
4. 将输出映射到 Network、lobe、broad anatomy 等可比层级。
5. 只报告 coarse consistency。

### 9.4 论文数字

当前 Table S4 submission route：

```text
n = 65
Network Top1 = 15.38%
Network Top3 = 40.00%
Broad anatomy Top1 = 13.85%
Broad anatomy Top3 = 64.62%
```

### 9.5 旧 MRI truth 数字是什么

曾出现过另一条 MRI truth / matching rule 路线：

```text
Network Top3 = 36.92%
Broad Top3 = 80.00%
```

这不是当前 Table S4 submission route。它只能列为 legacy/inconsistency 或 alternate matching rule，不能替换当前投稿数字。

## 10. GSE189919 projection feasibility

### 10.1 这条路线问什么问题

问题是：

> 外部 CSF RNA-seq 矩阵的 gene space 是否足够覆盖当前 projector/model genes，使得 projection 技术上可运行？

它不问：

> 预测脑区是否正确？

因为 GSE189919 没有 patient-level anatomical truth。

### 10.2 主要脚本

`scripts/analyze_gse189919_csf_tracing_validation.py`

### 10.3 算法步骤

1. 读取 GSE189919 count matrix 和软 metadata。
2. 标准化 sample/gene identifiers。
3. 计算 model/projector gene overlap。
4. 对样本生成 Network prediction distributions。
5. 输出 route agreement、broad category distribution、sample QC。
6. 可运行 optional algorithm audit。

### 10.4 论文数字

```text
Projector gene overlap = 15622 / 21668
Coverage = 72.10%
Conclusion = projection feasibility only
```

### 10.5 algorithm audit 失败如何解释

该脚本包含一个 optional legacy algorithm audit，用来检查一些旧 baseline route 的实现假设。如果它报：

```text
Algorithm audit failed:
baseline_production_tpm ...
baseline_production_count_cpm ...
```

这不等于投稿 accuracy route 失败。它说明旧 baseline/production comparison 的某些技术假设不再支持作为投稿准确率证据。

当前 workflow 应将其标为：

```text
EXPECTED_LEGACY_ALGORITHM_LIMITATION
```

而不是把 GSE189919 当成 localization accuracy。

## 11. Public audit 与 controlled-data validation

### 11.1 为什么需要双模式

Bo2023 raw expression matrix 是 controlled / non-public 数据。公开 GitHub release 不能直接包含完整 Bo2023 expression matrix。

因此 HTML workflow 设计为双模式：

1. Public audit mode
   - GitHub clone 后可运行。
   - 检查代码路径、source data、表格、测试、公开模型和 claim boundary。
   - 缺少 Bo2023 controlled files 时跳过 full region tracing test。

2. Controlled-data validation mode
   - 当本地存在授权 Bo2023/AHBA/TCGA/GSE 输入时自动启用。
   - 运行完整复现。

辅助脚本：

```text
scripts/stage_validation_data_under_tests.py
```

该脚本用于在 Victor 服务器上把授权验证数据集中复制到 `tests/controlled_data/` 下。它的边界是只允许从 `/storage/wangzhen` 允许范围内查找和复制文件，避免 workflow 为了复现而访问其他用户目录或不受控路径。

### 11.2 测试行为

公开 release 缺少 controlled Bo2023 raw files 时：

```text
16 passed, 1 skipped
```

Victor 或本地授权数据齐全时：

```text
17 passed
```

这两个结果都合理，取决于 controlled data 是否存在。

### 11.3 数据集中放置

Victor 上推荐将授权数据放在：

```text
/storage/wangzhen/cfrna-brain-tracing-0.1.6/tests/controlled_data/
```

其中 Bo2023 常用结构：

```text
tests/controlled_data/bo2023/
  mfas5_819samples_28415genes_featurecounts_counts.txt
  mfas5_819samples_23605genes_vsd4_rmbatch.xls
  Information of sequenced samples_update_full878_filter819.xlsx
  04_expressed_genes_neocortex_plus_subcortical.cleaned_symbols.csv
```

## 12. HTML workflow 如何复现

主 workflow：

```text
output/cfrna_braintrace_v016_jupyter_validation_workflow.html
```

它的作用不是发明新算法，而是把已锁定脚本按 Jupyter Lab cell 形式串起来：

1. 环境和路径检测。
2. public source data 检查。
3. `pytest tests`。
4. projected-VSD Network LOSO。
5. projected-VSD Network LOMO。
6. formal three-tier LOSO。
7. formal three-tier LOMO。
8. AHBA external validation。
9. TCGA/BraTS coarse consistency。
10. GSE189919 projection feasibility。
11. legacy/inconsistency route 标记。
12. 写出 final summary table。

每条路线应输出：

- command
- inputs
- outputs
- target numbers
- actual numbers
- consistency
- difference reason
- status
- result table
- result plot

最终汇总：

```text
output/v016_jupyter_validation/jupyter_validation_summary.csv
output/v016_jupyter_validation/jupyter_validation_summary.md
```

## 13. 怎样读 Top1、Top3 和 median true-rank

### 13.1 Top1

预测排名第一的标签等于 truth label。

例如：

```text
truth = Network A
prediction = [Network A, Network B, Network C]
Top1 hit = True
```

### 13.2 Top3

truth label 出现在前三个候选中。

例如：

```text
truth = Network A
prediction = [Network B, Network C, Network A]
Top1 hit = False
Top3 hit = True
```

Top3 更适合 candidate ranking system，因为工具的目的是给出候选来源，而不是宣称单一确定定位。

### 13.3 median true-rank

每个样本都有 truth label 在完整排序中的位置。取中位数就是 median true-rank。

如果 median true-rank = 1.0，表示至少一半样本的真实标签排在第一位。

## 14. 常见误解与正确解释

### 14.1 “Network n=819，region n=814/812，是不是漏掉样本？”

不是漏掉。原因是 region 层需要 truth region 在训练 reference 中存在。

如果 held-out 样本的 truth region 在训练折不存在，那么算法没有合法 reference 可以评价该 region。把它算错或算对都会制造假证据。

正确做法：

- Network 仍评价。
- group/exact 跳过，并在 unsupported samples 文件中披露。

### 14.2 “projected-VSD exact-region 低，是否推翻算法？”

不推翻。论文主线没有把 projected-VSD direct exact-region scoring 当作最终 route。

projected-VSD 的作用是 Network beam。exact-region 使用 logCPM-compatible local reranking，且是 exploratory ranking。

### 14.3 “GSE189919 有预测结果，为什么不算准确率？”

因为没有 anatomical truth。没有真值就不能算 accuracy。

它只能说明：

- gene overlap 是否足够；
- projection 是否能运行；
- 输出分布是否可审计。

### 14.4 “TCGA/BraTS 是 human glioma，为什么能用？”

它不是用来证明 macaque exact-region localization，而是 stress test / coarse consistency。

MRI-derived labels 与 Bo2023 macaque region label 不同。因此只报告 Network/broad anatomy 的 coarse consistency。

### 14.5 “旧数字还能不能引用？”

旧 baseline 或旧 denominator 数字不能作为当前 submission route。

允许出现的位置：

- audit changelog；
- legacy/inconsistency section；
- difference reason。

不应出现的位置：

- abstract 主结果；
- Table S2/S4 当前 route；
- Figure 1 当前 summary；
- app 当前 validation context。

## 15. 当前投稿数字索引

### 15.1 Bo2023 internal validation

| Dataset | Route | Endpoint | n | Top1 | Top3 |
|---|---|---:|---:|---:|---:|
| Bo2023 LOSO | Projected VSD Network | Network | 819 | 58.00% | 91.58% |
| Bo2023 LOMO | Projected VSD Network | Network | 819 | 53.72% | 91.33% |
| Bo2023 LOSO | Formal hybrid | Network | 819 | 58.24% | 92.19% |
| Bo2023 LOSO | Formal hybrid | Resolution group | 814 | 44.47% | 72.36% |
| Bo2023 LOSO | Formal hybrid | Exact region | 814 | 22.48% | 45.33% |
| Bo2023 LOMO | Formal hybrid | Network | 819 | 57.75% | 91.21% |
| Bo2023 LOMO | Formal hybrid | Resolution group | 812 | 41.38% | 69.09% |
| Bo2023 LOMO | Formal hybrid | Exact region | 812 | 22.17% | 42.36% |

### 15.2 External validation

| Dataset | Route | Endpoint | Top1 | Top3 | Claim boundary |
|---|---|---:|---:|---:|---|
| AHBA | Hybrid | Network | 74.68% | 94.42% | mapped-label transfer |
| AHBA | Hybrid | Resolution group | 36.26% | 67.03% | mapped-label transfer |
| AHBA | Hybrid | Exact region | 24.18% | 42.86% | exact-mapped subset |
| TCGA/BraTS | Hybrid | Network | 15.38% | 40.00% | coarse consistency |
| TCGA/BraTS | Hybrid | Broad anatomy | 13.85% | 64.62% | coarse consistency |
| GSE189919 | Projection feasibility | Gene overlap | 15622/21668 | 72.10% | no accuracy claim |

## 16. 审稿复现时应检查什么

人工审查建议按以下顺序：

1. 运行 `pytest tests`
   - public mode 期望 `16 passed, 1 skipped`
   - controlled-data mode 期望 `17 passed`

2. 检查 source tables
   - `TableS1_internal_validation_design.csv`
   - `TableS2_internal_validation_results.csv`
   - `TableS4_external_validation_results.csv`
   - `Figure1_validation_summary.csv`

3. 检查 Network 与 region 分母
   - LOSO Network: 819
   - LOSO group/exact: 814
   - LOMO Network: 819
   - LOMO group/exact: 812

4. 检查 unsupported samples 是否输出
   - LOSO: `hybrid_formal_loso_region_unsupported_samples.csv`
   - LOMO: `formal_lomo_region_unsupported_samples.csv`

5. 检查 legacy 数字是否只出现在 legacy/inconsistency 说明
   - `92.38%`
   - `36.92% / 80.00%`

6. 检查 GSE189919 是否没有 accuracy claim。

7. 检查 final summary 是否记录每条路线的：
   - command
   - inputs
   - outputs
   - target
   - actual
   - consistency
   - difference reason
   - status

## 17. 一句话理解整个项目

cfRNA-BrainTrace v0.1.6 不是声称从任意 cfRNA 样本精确定位到一个脑区；它是一个分层候选排序工具。当前投稿主线先用 projected-VSD 在粗 Network 层生成 Top3 候选束，再在该候选束内用 logCPM-compatible reference 做局部 group/exact ranking。验证结果显示 Network Top3 最稳健，resolution group 有中等支持，exact region 只能作为探索性候选排序。外部数据根据真值分辨率严格限制解释边界。
