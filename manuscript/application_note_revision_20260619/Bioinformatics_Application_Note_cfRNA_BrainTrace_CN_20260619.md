# cfRNA-BrainTrace：基于灵长类转录组图谱的 RNA-seq 分层脑来源推断

**文章类型：** Bioinformatics Application Note  
**栏目：** Gene expression  
**作者与单位：** [投稿前补充]  
**通讯作者：** [投稿前补充]

## 摘要

### Summary

cfRNA-BrainTrace 是一个 Python 与 Streamlit 软件，用于从 RNA 表达谱中进行分层脑来源候选排序。当前锁定路线将 query 从 logCPM/logTPM 兼容表达空间投影到 Bo2023-like VSD 空间，仅用于 10 类 Network Top3 候选生成；随后 resolution group 与 exact region 均回到 logCPM 兼容的局部表达空间重排序。该设计避免把投影误写成新 atlas，并把解剖分辨率作为显式输出。在 Bo2023 内部验证中，正式三级 hybrid 路线的 Network Top3 在 LOSO 和 LOMO 中分别为 0.9238 和 0.9121；resolution-group Top3 分别为 0.7236 和 0.6909；exact-region Top3 分别为 0.4533 和 0.4236。AHBA 外部映射标签验证支持 hybrid 路线，TCGA/BraTS 只支持粗粒度解剖一致性。该软件适用于可复现的分层候选排序、质量控制和溯源分辨率边界审计。

### Availability and Implementation

Python 3.11+；支持命令行和 Streamlit 网页界面。开发仓库目前为私有仓库 **https://github.com/wz7717/cfrna-brain-tracing**；投稿前需完成公开发布、归档 DOI 和许可证确认。稳定在线应用地址为 **[投稿前补充]**。

## 1 引言

cfRNA 或组织 RNA 表达谱中可能保留来源组织信息，但脑内解剖定位受到区域转录相似性、图谱分辨率和样本域偏移限制。若直接给出单一 exact-region 标签，容易把粗粒度可信信号误写成精确定位。cfRNA-BrainTrace 因此采用分层报告策略：Network Top3 是主 endpoint，resolution group 是主要 region-level 粒度，exact region 仅作为局部候选排序。

## 2 软件与方法

### 2.1 输入与预处理

软件接收含基因符号和表达量的表格输入，并保留可用的样本、受试者、诊断和解剖元数据。输入被转换为 logCPM/logTPM 兼容空间后进行基因匹配、覆盖度审计和范围检查。每个样本都会输出匹配基因数、非零覆盖度、entropy、margin 和 route 标识。

### 2.2 Reference projection 与三级推断

训练参考仍为原始 Bo2023 表达资源：Network beam 使用 Bo2023 VSD/batch-removed 参考，resolution group 与 exact region 使用 raw-feature-count-derived logCPM 参考。projector 只把被留出样本或外部 query 映射到 Bo2023-like VSD 空间，用于 Network Top3 beam。projector 的 sample-level 拟合很强（global Pearson 0.9961，median sample Pearson 0.9969），但 gene-level 拟合较弱，因此不能写成“每个基因都被精确恢复”。

Network Top3 beam 生成后，候选脑区被限制在这些 Network 内；resolution group 与 exact region 再使用 logCPM 局部判别基因重排序。这个 hybrid 设计来自直接 exact-region 诊断结果：projected VSD 在 exact 层并非始终优于其他路线，尤其在 LOMO 中低于 native VSD。

## 3 实现

主线代码已经接入正式三级路线：Network projector route 位于 `core/network_tracing.py`，three-tier region route 位于 `core/bo2023_region_tracing.py`，网页触发和 route 说明位于 `app/pages/tracing_page.py`。数据库样本 `19R348` smoke test 确认输出 route 为 `projected_vsd_network_top3_logcpm_resolution_local_exact`，Network overlap 为 199/200，Region reference source 为 raw featurecounts logCPM。

## 4 验证结果

### 4.1 Network 层内部验证

在 819 个 Bo2023 样本的 fold-local LOSO Network 验证中，projected VSD 的 Network Top1/Top3 为 0.5800/0.9158，高于 logCPM baseline 和 native VSD。在 strict LOMO 验证中，projected VSD 的 Network Top1/Top3 为 0.5372/0.9133，同样高于两个对照路线。这是把 projected VSD 用作 Network candidate beam 的主要证据。

### 4.2 Direct exact-region 仅作为诊断基线

direct projected-VSD exact scoring 在 LOSO 中的 exact Top3 为 0.3563，高于 logCPM baseline；但在 LOMO 中 exact Top3 为 0.3116，低于 native VSD 的 0.3473。因此不能写“projected VSD 在 exact-region 层全面优于其他路线”，direct exact 也不应作为主线 endpoint。

### 4.3 正式三级 hybrid 验证

完整三级路线为 projected VSD Network Top3 beam -> logCPM resolution group rerank -> logCPM local exact rerank。LOSO 中 Network、resolution group、exact region 的 Top3 分别为 0.9238、0.7236、0.4533；LOMO 中分别为 0.9121、0.6909、0.4236。Group Top3 明显高于 Exact Top3，说明 resolution group 是更稳健的 region-level 报告粒度。

### 4.4 外部验证

AHBA human brain RNA-seq 映射标签验证中，hybrid 的 Network Top1/Top3 为 0.7468/0.9442，Group Top1/Top3 为 0.3626/0.6703，Exact Top1/Top3 为 0.2418/0.4286。hybrid 的 Exact Top3 高于 logCPM baseline 和 projected-VSD-only，支持“projected VSD 保留 Network 优势，logCPM downstream 修复 group/exact 层损失”的解释。

TCGA/BraTS glioma tissue RNA-seq + MRI label 只能用于 coarse anatomical consistency，因为 MRI 真值是人脑 atlas 标签，不是 Bo2023 macaque exact-region ID。hybrid 的 Network Top3 为 0.4000，Broad anatomy Top3 为 0.6462，均高于 logCPM baseline。GSE189919 只证明外部矩阵可投影到 projector gene space，不能写成 brain-origin accuracy 验证。

## 5 结论

cfRNA-BrainTrace 将三层 atlas-based RNA 溯源路线封装为可复现的软件。当前证据最支持 Network Top3 和 resolution-group 层级报告；exact-region Top1 仍受图谱分辨率和相邻脑区转录相似性限制。软件应用边界应定位为分层候选排序、质量控制和 claim-boundary 审计，而不是无标签 biofluid 的临床精确定位。

## 图注

**图1. cfRNA-BrainTrace 三级路线与验证。** query 仅在 Network beam 阶段投影到 Bo2023-like VSD；resolution group 和 exact region 回到 logCPM 局部空间重排序。内部 LOSO/LOMO 支持 Network Top3 与 group-level endpoint；AHBA 支持跨物种映射标签验证；TCGA/BraTS 仅支持粗粒度解剖一致性。
