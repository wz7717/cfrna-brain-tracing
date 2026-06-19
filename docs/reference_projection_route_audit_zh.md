# Reference-projected VSD 路线审计

## 审计范围

本审计覆盖本线程新增的 reference-projected VSD 分支、gene symbol 清理、内部 Bo2023 验证、AHBA 外部验证、TCGA/BraTS MRI 标注外部验证，以及 PDF 汇总报告生成流程。

## 推荐主线

当前推荐的跨域主线是：

`projected VSD Network Top3 beam -> logCPM resolution group -> logCPM local exact rerank`

理由：

- `projected VSD` 在跨域输入上更适合作为 Network Top3 候选入口。
- `logCPM` 在 region-local rerank 中更稳，能修复 pure projected VSD 在 exact region 层的损失。
- 完整三级 LOSO 中，该 hybrid 路线达到 `Network Top3=0.9238`、`Group Top3=0.7236`、`Exact Top3=0.4533`。
- AHBA 外部正式三级中，hybrid 达到 `Exact Top3=0.4286`，明显高于 pure projected VSD 和 logCPM baseline。

## 关键结果

| 验证 | Network Top3 | Group Top3 | Exact Top1 | Exact Top3 | 结论 |
| --- | ---: | ---: | ---: | ---: | --- |
| 内部完整三级 LOSO hybrid | 0.9238 | 0.7236 | 0.2248 | 0.4533 | 当前内部最支持 hybrid 的结果 |
| 内部完整三级 LOMO hybrid | 0.9121 | 0.6909 | 0.2217 | 0.4236 | exact Top3 与 logCPM 持平，修复 pure projected 损失 |
| AHBA formal hybrid | 0.9442 | 0.6703 | 0.2418 | 0.4286 | 外部 mapped exact 上最优 |
| TCGA/BraTS labeled hybrid | 0.4000 | NA | NA | NA | MRI truth 是 human atlas 标签，只报告 network/lobe/broad |

## 非推荐路线

- `pure projected VSD exact rerank`：内部 LOMO exact Top3 为 `0.3793`，低于 hybrid `0.4236`。
- `direct exact region scoring`：未经过 formal network beam 和 local rerank，exact 层偏弱，只保留为诊断基线。
- `TCGA/BraTS Bo2023 exact accuracy`：MRI truth 使用 human atlas label，不是 Bo2023 macaque region ID，不能报告严格 Bo2023 exact accuracy。

## 必要修复

- `core/reference_projection.py` 提供 `clean_excel_date_gene_symbol()`，用于清理 Excel 自动日期转换造成的 gene symbol 污染。
- `scripts/run_bo2023_v2_loso_validation.py` 已接入该清理函数，避免主线 LOSO 映射再次引入日期样式 gene symbol。
- 新增 gene symbol 审计与清理脚本用于复跑和追踪变更。

## 入仓策略

应入仓：

- 可复跑代码：`core/reference_projection.py` 与新增 `scripts/*.py`。
- 路线计划与审计文档。
- 最终 PDF 报告和报告生成脚本。

不应入仓：

- `results/` 下的大型中间结果，已由 `.gitignore` 忽略。
- PDF 渲染检查 PNG，属于临时 QA 产物。
- 未经确认的生产模型重生成产物，除非明确作为 cleaned production model 发布。

## 后续建议

1. 将 hybrid 路线作为外部验证默认路径。
2. 对 AHBA 的 Putamen/Caudate/Insula 等标签继续做误差审计。
3. 建立 human atlas label 到 Bo2023 region group 的正式映射层，再扩展 TCGA/BraTS exact 层评估。
4. 若要替换生产模型，应单独开一次模型发布流程，固定训练输入、生成命令、指标和校验摘要。
