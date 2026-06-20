# cfRNA-BrainTrace 中文补充材料

## 补充方法

### S1 正式验证路线

所有验证均按照软件正式路线组织。query profile 先表示为 logCPM/logTPM 兼容表达空间，再仅在 Network scoring 阶段投影到 Bo2023-like VSD 空间。Top3 Networks 构成 candidate beam。随后 resolution-group 和 exact-region candidates 都在该 beam 内使用 logCPM 兼容局部表达空间重排序。该设计把 broad candidate generation 与 fine regional interpretation 分开。

### S2 内部验证设计

内部验证使用 Bo2023 猕猴脑参考。Leave-one-sample-out validation 测试在留出单个样本时能否恢复其标签；leave-one-monkey-out validation 每次留出一只猴的全部样本，用于测试 donor-level generalization。指标分别在 Network、resolution group 和 exact region 三个层级报告。Top1、Top3 和 median true-rank 用于区分单一标签命中与候选列表召回。

### S3 外部验证设计

外部验证按数据集标签分辨率解释。AHBA human brain RNA-seq 用于 mapped-label validation，因为其解剖标签可映射到 Network、resolution group 和部分 exact-region 标签。TCGA/BraTS glioma tissue RNA-seq + MRI-derived labels 只用于 coarse anatomical consistency，因为其真值是人脑影像标签，不是 Bo2023 猕猴 exact-region ID。GSE189919 用于测试外部矩阵是否可投影到模型基因空间；由于缺少 patient-level anatomical truth，不用于 accuracy estimation。

### S4 图表组织

补充图 S1 给出主路线与验证总结图。补充表 S1-S5 分别给出内部验证设计、内部验证结果、外部验证设计、外部验证结果以及图表索引。这些材料共同说明正式路线如何验证、每个验证设置可支持何种结论。

## 补充结果

### S1 内部路线选择

第一组内部分析测试 Network-level candidate generation。projected VSD 在 LOSO 中 Network Top1/Top3 为 58.00%/91.58%，在 LOMO 中为 53.72%/91.33%，在 Network Top3 上均高于 logCPM baseline 与 native VSD。Direct exact-region scoring 的结果更低且更不稳定，尤其在 LOMO 中表现不足。因此正式路线使用 projected VSD 生成 broad Network beam，并使用 logCPM 兼容表达进行下游 local reranking。

### S2 正式内部三级验证

完整 LOSO validation 中，正式路线的 Network、resolution-group 和 exact-region Top3 分别为 92.38%、72.36% 和 45.33%。完整 LOMO validation 中，对应 Top3 分别为 91.21%、69.09% 和 42.36%。Median true-rank 从 Network 到 exact-region 层级逐步增大，符合更细解剖分辨率下不确定性升高的预期。因此 resolution group 是更适合的 region-level endpoint，exact-region output 保留为候选排序。

### S3 外部验证结果

AHBA 中，hybrid route 的 Network Top1/Top3 为 74.68%/94.42%，resolution-group Top1/Top3 为 36.26%/67.03%，exact-region Top1/Top3 为 24.18%/42.86%。Hybrid exact Top3 高于 logCPM baseline 的 30.77% 和 projected-VSD-only 的 29.67%。TCGA/BraTS 中 hybrid Network Top3 为 40.00%，broad-anatomy Top3 为 64.62%，只支持 coarse anatomical consistency。GSE189919 覆盖 15,622/21,668 projector genes，即 72.10% gene-space coverage，支持 projection feasibility，而不是 source-localization accuracy。
