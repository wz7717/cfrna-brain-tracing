# Jupyter Lab 人工验证 workflow 使用手册

本文档专门说明如何使用：

```text
output/cfrna_braintrace_v016_jupyter_validation_workflow.html
```

在 Jupyter Lab 中逐步复现 cfRNA-BrainTrace v0.1.6 的人工验证流程。

如果只记一句话：

```text
这个 HTML workflow 是审稿复现入口；它按 public audit / controlled-data validation 双模式运行，每个验证 block 都会生成命令记录、结果表、结果图和最终 summary，并严格区分当前投稿路线与 legacy/inconsistency。
```

## 1. workflow 文件在哪里

主文件：

```text
output/cfrna_braintrace_v016_jupyter_validation_workflow.html
```

建议在浏览器中打开这个 HTML，然后把每个 code block 复制到 Jupyter Lab notebook 中顺序运行。

它不是新的算法实现，而是把仓库里已经锁定的脚本串起来：

```text
pytest tests
Bo2023 projected-VSD Network LOSO
Bo2023 projected-VSD Network LOMO
Bo2023 formal three-tier LOSO
Bo2023 formal three-tier LOMO
AHBA mapped-label external validation
TCGA/BraTS coarse consistency
GSE189919 projection feasibility
GSE optional technical audit
final summary export
```

## 2. 运行前准备

### 2.1 项目根目录

Victor 服务器推荐路径：

```text
/storage/wangzhen/cfrna-brain-tracing-0.1.6
```

本地 Windows 当前仓库路径：

```text
D:\Download\cfrna-brain-tracing-streamlit-cloud-ready
```

workflow 会优先读环境变量：

```text
CFRNA_BRAINTRACE_ROOT
```

如果在 Victor 上运行，建议先设置：

```bash
export CFRNA_BRAINTRACE_ROOT=/storage/wangzhen/cfrna-brain-tracing-0.1.6
```

或在 notebook 里确认 `ROOT` 自动识别到了正确目录。

### 2.2 Python 环境

HTML 启动页列出需要的软件包：

```text
numpy
pandas
scipy
matplotlib
openpyxl
nibabel
pytest
jupyterlab
nbformat
nbclient
```

如果缺包，可安装：

```bash
python -m pip install numpy pandas scipy matplotlib openpyxl nibabel pytest jupyterlab nbformat nbclient
```

Victor 上如果已有环境：

```text
/home/wangzhen/.conda/envs/cfrna_bt_v016/bin/python
```

优先使用该环境。

## 3. public audit 与 controlled-data validation

workflow 有双模式。

### 3.1 Public audit mode

公开 GitHub clone 后可运行。

特点：

- 不要求 Bo2023 controlled raw expression matrix。
- 可以检查代码路径、公开模型、论文 source data、Table S1-S6、Figure source data。
- 遇到需要 controlled data 的 full validation block 会 skip。

pytest 预期：

```text
16 passed, 1 skipped
```

这个结果是合理的。

### 3.2 Controlled-data validation mode

如果本地存在授权数据，workflow 会自动进入 controlled-data validation。

特点：

- 运行 Bo2023、AHBA、TCGA/BraTS、GSE189919 的完整复现 block。
- 输出每条路线的实际数字和图表。

pytest 预期：

```text
17 passed
```

### 3.3 controlled data 推荐位置

workflow 默认检查：

```text
tests/controlled_data/
```

典型结构：

```text
tests/controlled_data/bo2023/
tests/controlled_data/ahba/
tests/controlled_data/tcga_brats/
tests/controlled_data/gse189919/
```

Bo2023 关键文件：

```text
tests/controlled_data/bo2023/mfas5_819samples_28415genes_featurecounts_counts.txt
tests/controlled_data/bo2023/mfas5_819samples_23605genes_vsd4_rmbatch.xls
tests/controlled_data/bo2023/Information of sequenced samples_update_full878_filter819.xlsx
tests/controlled_data/bo2023/04_expressed_genes_neocortex_plus_subcortical.cleaned_symbols.csv
```

## 4. 输出目录

workflow 的主输出目录：

```text
output/v016_jupyter_validation/
```

每个 route 会在下面创建自己的子目录，例如：

```text
output/v016_jupyter_validation/pytest_tests/
output/v016_jupyter_validation/bo2023_projected_vsd_network_loso/
output/v016_jupyter_validation/bo2023_projected_vsd_network_lomo/
output/v016_jupyter_validation/bo2023_formal_three_tier_loso/
output/v016_jupyter_validation/bo2023_formal_three_tier_lomo/
output/v016_jupyter_validation/ahba_mapped_label_external_validation/
output/v016_jupyter_validation/tcga_brats_coarse_consistency/
output/v016_jupyter_validation/gse189919_projection_feasibility/
```

最终汇总：

```text
output/v016_jupyter_validation/jupyter_validation_summary.csv
output/v016_jupyter_validation/jupyter_validation_summary.md
```

每个 block 通常还会输出：

```text
*_display_table.csv
*.png
run.log
```

## 5. Cell 1：路径、帮助函数和输出目录

这个 cell 做几件事：

1. 识别 `ROOT`。
2. 把 `ROOT` 加入 `sys.path`。
3. 创建 `AUDIT_DIR`。
4. 定义 `run_cmd`、`missing_paths`、`skip_controlled` 等帮助函数。

你要检查：

```text
ROOT 是否是当前项目根目录
AUDIT_DIR 是否在 output/v016_jupyter_validation
```

如果 `ROOT` 指错了，后面所有输入路径都会错。

## 6. Cell 1b：展示表和图形工具

这个 cell 定义：

```text
save_display_table
save_metric_barplot
metric_table
pct
```

作用：

- 把每条验证路线的结果保存成展示表。
- 把 Top1/Top3 或 gene overlap 保存成 PNG 图。

如果这个 cell 没跑，后面会出现：

```text
NameError: save_display_table is not defined
```

## 7. Cell 2：软件包和 Jupyter 环境检查

这个 cell 检查 required packages。

如果缺包，先安装再继续。

常见缺包：

```text
openpyxl
nibabel
pytest
```

`openpyxl` 用于读 Excel metadata。

## 8. Cell 3：论文目标数字和 source data 读取

这个 cell 把当前投稿目标数字写进 notebook。

目标数字包括：

```text
Bo2023 projected-VSD Network LOSO 58.00% / 91.58%
Bo2023 projected-VSD Network LOMO 53.72% / 91.33%
Formal LOSO Network 58.24% / 92.19%
Formal LOSO group/exact 72.36% / 45.33%
Formal LOMO Network 57.75% / 91.21%
Formal LOMO group/exact 69.09% / 42.36%
AHBA Network 74.68% / 94.42%
TCGA/BraTS Network Top3 40.00%
TCGA/BraTS broad Top3 64.62%
GSE189919 overlap 15622 / 21668, 72.10%
```

这个 cell 的作用是让后面每条路线都有 target 可比。

## 9. Cell 4：双模式输入文件检查

这个 cell 定义：

```text
CONTROLLED_DATA_DIR = ROOT / "tests" / "controlled_data"
INPUTS
CONTROLLED_INPUTS
DATA_AVAILABLE
```

它会判断哪些输入存在。

如果缺少 Bo2023 controlled files，你会看到类似：

```text
SKIP controlled-data validation: Bo2023 projected-VSD Network LOSO
missing: ...
```

这不是 workflow 错误，而是 public audit mode 的正常行为。

## 10. Cell 5：运行 pytest tests

命令：

```text
python -m pytest tests -p no:cacheprovider --basetemp ...
```

### 10.1 public mode 预期

```text
16 passed, 1 skipped
```

原因：

```text
公开 release 缺少 controlled Bo2023 raw expression inputs，full region tracing test skip。
```

### 10.2 controlled mode 预期

```text
17 passed
```

原因：

```text
tests/controlled_data/bo2023/ 中存在授权 Bo2023 输入，full region tracing test 可以运行。
```

### 10.3 常见 AssertionError

如果 pytest 输出是 `17 passed`，但 cell 进入 public branch 断言 `16 passed, 1 skipped`，通常说明：

```text
cell 判断 controlled input 的路径与 pytest 实际使用的路径不一致。
```

解决思路：

- 检查 `INPUTS["Bo2023 counts"]`
- 检查 `INPUTS["Bo2023 sample info"]`
- 检查 gene map fallback 路径
- 确认 `ROOT` 是否正确

## 11. Cell 6：projected-VSD 使用边界检查

这个 cell 检查代码主线是否符合论文：

```text
projected-VSD 只用于 Network Top3 beam
resolution-group / exact-region 使用 logCPM-compatible local reranking
```

你要确认它没有把旧 baseline 或 direct projected-VSD exact route 当作当前 submission route。

## 12. Bo2023 projected-VSD Network LOSO block

HTML 章节：

```text
Bo2023 projected-VSD Network LOSO
```

运行脚本：

```text
scripts/run_bo2023_projected_vsd_loso.py
```

目标数字：

```text
n = 819
Top1 = 58.00%
Top3 = 91.58%
```

输出：

```text
bo2023_projected_vsd_loso_detail.csv
bo2023_projected_vsd_loso_route_summary.csv
bo2023_projected_vsd_loso_summary.json
bo2023_projected_vsd_loso_display_table.csv
bo2023_projected_vsd_loso_metrics.png
run.log
```

解释：

```text
这是 Network-only 方法合理性验证，不是完整 three-tier route。
```

## 13. Bo2023 projected-VSD Network LOMO block

运行脚本：

```text
scripts/run_bo2023_projected_vsd_lomo.py
```

目标数字：

```text
n = 819
Top1 = 53.72%
Top3 = 91.33%
```

输出：

```text
bo2023_projected_vsd_lomo_detail.csv
bo2023_projected_vsd_lomo_folds.csv
bo2023_projected_vsd_lomo_route_summary.csv
bo2023_projected_vsd_lomo_summary.json
bo2023_projected_vsd_lomo_display_table.csv
bo2023_projected_vsd_lomo_metrics.png
```

解释：

```text
这是独立 projected-VSD Network LOMO，不等同于 formal LOMO Network。
```

## 14. Bo2023 formal three-tier LOSO block

运行脚本：

```text
scripts/run_bo2023_hybrid_formal_loso.py
```

目标数字：

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

关键输出：

```text
hybrid_formal_loso_network_route_metrics.csv
hybrid_formal_loso_resolution_group_route_metrics.csv
hybrid_formal_loso_exact_region_route_metrics.csv
hybrid_formal_loso_region_unsupported_samples.csv
hybrid_formal_loso_summary.json
hybrid_formal_loso_display_table.csv
hybrid_formal_loso_metrics.png
```

重点检查：

```text
Network n = 819
group/exact n = 814
unsupported = 5
```

不要把旧 `92.38%` 当当前 LOSO Network Top3。

## 15. Bo2023 formal three-tier LOMO block

运行脚本：

```text
scripts/run_bo2023_projected_vsd_formal_lomo.py
```

目标数字：

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

关键输出：

```text
formal_lomo_network_route_metrics.csv
formal_lomo_resolution_group_route_metrics.csv
formal_lomo_exact_region_route_metrics.csv
formal_lomo_region_unsupported_samples.csv
formal_lomo_validation_summary.json
formal_lomo_display_table.csv
formal_lomo_metrics.png
```

重点检查：

```text
Network n = 819
group/exact n = 812
unsupported = 7
```

## 16. AHBA mapped-label external validation block

运行脚本：

```text
scripts/run_ahba_projected_vsd_formal_three_tier_external.py
```

目标数字：

```text
Network Top1/Top3 = 74.68% / 94.42%
Resolution group Top1/Top3 = 36.26% / 67.03%
Exact region Top1/Top3 = 24.18% / 42.86%
```

关键输出：

```text
ahba_formal_three_tier_sample_detail.csv
ahba_formal_three_tier_metrics.csv
ahba_formal_three_tier_special_labels.csv
ahba_formal_three_tier_summary.json
ahba_formal_three_tier_display_table.csv
ahba_formal_three_tier_metrics_display.png
```

解释边界：

```text
AHBA = mapped-label external validation
不是 direct human-macaque exact anatomical equivalence。
```

## 17. TCGA/BraTS coarse consistency block

运行脚本：

```text
scripts/run_tcga_labeled_hybrid_formal_external.py
```

目标数字：

```text
Network Top1/Top3 = 15.38% / 40.00%
Broad anatomy Top1/Top3 = 13.85% / 64.62%
```

关键输出：

```text
tcga_labeled_hybrid_formal_sample_detail.csv
tcga_labeled_hybrid_formal_metrics.csv
tcga_labeled_hybrid_formal_summary.json
tcga_labeled_hybrid_formal_display_table.csv
tcga_labeled_hybrid_formal_metrics_display.png
```

解释边界：

```text
TCGA/BraTS = coarse anatomical consistency only
不是 Bo2023 exact-region localization。
```

旧 side route：

```text
Network Top3 = 36.92%
Broad Top3 = 80.00%
```

只能作为 legacy/inconsistency，不能替换当前 Table S4。

## 18. GSE189919 projection feasibility block

这个 block 复现 Table S4 gene-space overlap。

目标数字：

```text
Projector gene overlap = 15622 / 21668
Coverage = 72.10%
```

关键输出：

```text
gse189919_projection_feasibility_display_table.csv
gse189919_projection_feasibility_metrics.png
```

解释边界：

```text
GSE189919 没有 patient-level anatomical truth。
因此只支持 projection feasibility，不支持 localization accuracy。
```

## 19. Optional GSE 技术审计 block

运行脚本：

```text
scripts/analyze_gse189919_csf_tracing_validation.py
```

这个 block 是可选技术审计，不是投稿准确率路线。

如果出现：

```text
Algorithm audit failed:
baseline_production_tpm ...
baseline_production_count_cpm ...
```

workflow 应记录为：

```text
EXPECTED_LEGACY_ALGORITHM_LIMITATION
```

正确解释：

```text
旧 baseline/production comparison 的某些技术假设不应作为当前 submission accuracy evidence。
```

不要解释成：

```text
GSE189919 localization accuracy failed。
```

## 20. Final summary block

最后一个 cell 写出：

```text
jupyter_validation_summary.csv
jupyter_validation_summary.md
```

summary 每行应包含：

```text
route
status
command
inputs
output
target
actual
consistent
difference_reason
```

这是最终审稿复现报告的核心文件。

## 21. status 怎么读

### 21.1 PASS

表示该 route 已运行，实际数字与目标数字一致。

### 21.2 SKIPPED_CONTROLLED_DATA

表示缺少本地授权数据，public audit mode 下跳过。

这不是失败。

### 21.3 EXPECTED_LEGACY_ALGORITHM_LIMITATION

表示 optional GSE legacy technical audit 触发了已知旧路线限制。

这不是当前投稿准确率路线失败。

### 21.4 LEGACY_INCONSISTENCY

表示检测到旧数字或旧 source-data route，只能作为 legacy/inconsistency 说明。

不能引用为当前 submission route。

### 21.5 NOT_PRESENT_IN_CURRENT_RELEASE

表示旧 source-data 目录在当前 release 中不存在。

这通常是合理状态，因为当前 release 不应把旧 baseline 当主路线。

## 22. 常见错误和解释

### 22.1 `KeyError: 'Bo2023 gene map'`

原因：

```text
workflow 中某个 cell 用 INPUTS["Bo2023 gene map"]，但 INPUTS 没有这个 key。
```

修正后的逻辑应使用：

```text
INPUTS.get("Bo2023 gene map", fallback_path)
```

### 22.2 `SyntaxError: unterminated string literal`

原因通常是手工复制 cell 时重复粘贴或引号断裂。

处理：

```text
重新从 HTML 复制完整 cell，不要只复制半行。
```

### 22.3 pytest 是 `17 passed` 但断言 public skip

原因：

```text
controlled data 实际存在，但 workflow 的 missing_paths 检查走了 public branch。
```

处理：

```text
检查 ROOT
检查 CONTROLLED_DATA_DIR
检查 Bo2023 counts/sample-info/gene-map 路径
```

### 22.4 GSE algorithm audit failed

正确状态：

```text
EXPECTED_LEGACY_ALGORITHM_LIMITATION
```

不是当前投稿路线失败。

### 22.5 Bo2023 block 显示 skip

如果你希望完整复现而不是 public audit，说明 controlled files 没放到 workflow 能识别的位置。

检查：

```text
tests/controlled_data/bo2023/
```

## 23. 每次跑完后要保存什么

建议保存或提交给审稿人的文件：

```text
output/v016_jupyter_validation/jupyter_validation_summary.csv
output/v016_jupyter_validation/jupyter_validation_summary.md
各 route 的 *_display_table.csv
各 route 的 *.png
各 route 的 run.log
unsupported samples CSV
```

尤其是：

```text
hybrid_formal_loso_region_unsupported_samples.csv
formal_lomo_region_unsupported_samples.csv
```

它们证明 region 分母处理是透明的。

## 24. 最短运行顺序

如果你只想确认当前投稿版是否可复现，按这个顺序：

```text
1. 设置 CFRNA_BRAINTRACE_ROOT
2. 打开 HTML workflow
3. 复制 Cell 1 到 Jupyter Lab 并运行
4. 运行 Cell 1b
5. 运行 Cell 2
6. 运行 Cell 3
7. 运行 Cell 4
8. 运行 Cell 5 pytest
9. 顺序运行每个 validation block
10. 运行 final summary block
11. 打开 jupyter_validation_summary.csv 检查 status / target / actual
```

## 25. 最短复述

这个 HTML workflow 的使用原则是：

```text
先确认 ROOT 和环境。
再确认 public/controlled 模式。
pytest 通过后，逐条运行验证路线。
每条路线都要看 display table、PNG、run.log 和 target/actual 是否一致。
缺 controlled data 的 skip 是合理 public audit。
GSE optional audit failed 是 legacy limitation，不是 accuracy failure。
legacy 数字只能出现在 inconsistency 说明中。
最终以 jupyter_validation_summary.csv/.md 作为人工验证报告。
```
