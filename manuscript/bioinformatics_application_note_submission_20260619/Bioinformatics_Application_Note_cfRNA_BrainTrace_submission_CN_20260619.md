# cfRNA-BrainTrace：基于灵长类转录组图谱的 RNA-seq 分层脑来源推断

**文章类型：** Application Note  
**栏目：** Gene expression  
**作者与单位：** 作者信息将在投稿系统中提供  
**通讯作者：** 通讯作者信息将在投稿系统中提供

## 摘要

### Summary

cfRNA-BrainTrace 是一个 Python 与 Streamlit 应用，用于从 RNA 表达谱中进行分层脑来源候选排序。正式路线只在 Network Top3 beam 生成阶段使用 projected-VSD 表达，随后回到 logCPM 兼容表达空间进行 resolution-group 和 exact-region 重排序。软件同时报告解剖分辨率、置信度和适用边界，避免把稳定的粗粒度信号误写成确定性 exact-region 定位。

### Availability and Implementation

软件基于 Python 3.11+ 实现，支持命令行和 Streamlit 网页界面。源码托管于开发仓库 **https://github.com/wz7717/cfrna-brain-tracing**。最终投稿包将配套公开版本、归档 DOI、许可证声明和在线服务地址。

### Contact

通讯作者联系方式将在投稿系统中提供。

### Supplementary information

补充材料在线提供。

## 正文

### 研究动机

RNA 表达谱可能保留组织来源信息，但脑内溯源受到相邻脑区转录相似性、图谱粒度以及组织、肿瘤和 biofluid 样本域偏移的限制。对于 atlas-based brain RNA tracing，一个粗粒度候选可能较稳定，但 exact-region localization 未必成立。如果只输出单一 exact-region 标签，容易暗示数据并不支持的解剖分辨率。cfRNA-BrainTrace 因此返回 Network、resolution group 和 exact region 三个层级的候选结果，并同步给出置信度、覆盖度和适用范围诊断。

### 系统与方法

正式路线来自三组验证现象。第一，projected-VSD query representation 对 broad Network candidate generation 最稳定。第二，直接在 projected-VSD 空间做 exact-region scoring 并不总是优于 native VSD 或 logCPM-based alternatives。第三，在 Network beam 内进行 local reranking 能提高下游解剖层级的解释性。因此最终路线是一个 hybrid procedure，而不是单一表达尺度上的 classifier。

在实际操作中，cfRNA-BrainTrace 先将上传表达表标准化为 logCPM/logTPM 兼容空间，并与参考基因面板对齐。随后 query 被投影到 Bo2023-like VSD 空间，仅用于 10-class Network scoring。Top3 Networks 构成 candidate beam，beam 之外的候选脑区不进入下游 regional ranking。在保留的 beam 内，resolution group 使用 logCPM 兼容的局部表达空间重排序，exact region 再作为较低置信度的局部候选排序输出。该路线使用 Bo2023 猕猴脑转录组参考，但不把整个 atlas 转换成新的 projected atlas。主线实现位于 `core/network_tracing.py` 和 `core/bo2023_region_tracing.py`，命令行和 Streamlit 调用同一评分核心。

### 验证与输出

验证设计与软件路线保持一致。首先在 Network 层测试 projected-VSD query 是否适合作为 broad candidate beam。在 fold-local leave-one-sample-out Network 验证中，projected VSD 的 Top1/Top3 为 58.00%/91.58%，高于 logCPM baseline 和 native VSD。在 strict leave-one-monkey-out 验证中，projected VSD 达到 53.72%/91.33%。相比之下，direct exact-region projected-VSD scoring 并非在所有设置中最优，尤其在 LOMO 中不如 native VSD，因此 exact-region scoring 没有被设置为唯一 endpoint。

随后对完整正式三级路线进行端到端验证：projected-VSD Network Top3 beam generation，logCPM-compatible resolution-group reranking，以及 logCPM-compatible exact-region reranking。该路线在 LOSO 中的 Network Top3 为 92.38%，在 LOMO 中为 91.21%。Resolution-group Top3 分别为 72.36% 和 69.09%，而 exact-region Top3 分别为 45.33% 和 42.36%。这种随解剖分辨率升高而下降的准确率梯度支持把 Network Top3 作为主 endpoint，把 resolution group 作为更稳健的 region-level 输出；exact-region 结果保留为探索性局部候选排序，而不是确定性定位。

外部分析按照标签能支持的分辨率解释。Allen Human Brain Atlas 映射标签验证中，hybrid 路线的 Network Top1/Top3 为 74.68%/94.42%，resolution-group Top1/Top3 为 36.26%/67.03%，exact Top1/Top3 为 24.18%/42.86%；hybrid exact Top3 高于 logCPM baseline 和 projected-VSD-only，支持“projected VSD 用于 broad beam、logCPM 用于 local reranking”的路线选择。TCGA/BraTS glioma tissue RNA-seq + MRI label 只支持 coarse anatomical consistency，因为 MRI 真值是人脑 atlas 标签，而不是 Bo2023 猕猴 exact-region ID。在该数据集中，hybrid Network Top3 为 40.00%，broad-anatomy Top3 为 64.62%。GSE189919 用于验证 projection feasibility，而不是 accuracy，因为缺少可计算 brain-origin accuracy 的解剖真值。

### 使用边界

软件导出 ranked candidates、matched-gene coverage、entropy、margin、route identifier 和 resolution-specific warning。只有在存在独立解剖真值时才计算 accuracy。无标签 biofluid 预测作为 transfer stress test 或假设生成候选排序报告，而不是临床定位。最终 release package 计划包含非商业访问、测试数据、稳定 URL、版本化源码和归档 DOI。

## Funding

基金信息将在投稿系统中声明。

## Conflict of Interest

无。

## Data availability

Bo2023 猕猴脑转录组图谱、Allen Human Brain Atlas、TCGA/BraTS 和 GEO 数据集均可从原始仓库获取，并受其原始访问条件约束。可公开的 processed evaluation tables 和 figure source data 已按带版本的公开 release 与归档准备。

## References

Bakas,S. *et al.* (2017) Advancing The Cancer Genome Atlas glioma MRI collections with expert segmentation labels and radiomic features. *Sci. Data*, **4**, 170117.

Bo,T. *et al.* (2023) Brain-wide and cell-specific transcriptomic insights into MRI-derived cortical morphology in macaque monkeys. *Nat. Commun.*, **14**, 1283.

Hawrylycz,M.J. *et al.* (2012) An anatomically comprehensive atlas of the adult human brain transcriptome. *Nature*, **489**, 391-399.

Vorperian,S.K. *et al.* (2022) Cell types of origin of the cell-free transcriptome. *Nat. Biotechnol.*, **40**, 855-861.

## 图注

**图 1. cfRNA-BrainTrace 三级路线与验证证据。** query 只在 Network Top3 beam 阶段投影到 Bo2023-like VSD 空间；下游 resolution-group 与 exact-region reranking 在 logCPM 兼容的局部表达空间中完成。内部验证支持 LOSO 和 LOMO 设置下较高的 Network Top3 accuracy，同时 resolution-group accuracy 持续高于 exact-region accuracy。AHBA 支持 mapped-label 外部验证，而 TCGA/BraTS 只支持 coarse anatomical consistency。

**Alt text：** 多面板图概述 cfRNA-BrainTrace 工作流与验证结果。流程图显示 query 预处理、projected-VSD Network Top3 beam 生成和 logCPM 局部重排序。柱状图显示内部验证 Network Top3 接近 92%，resolution-group Top3 约 69%-72%，exact-region Top3 约 42%-45%。外部面板显示 AHBA Network Top3 为 94.42%，TCGA/BraTS broad-anatomy Top3 为 64.62%。
