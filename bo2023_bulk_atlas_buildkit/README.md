# Bo2023 supplementary build kit

这套文件来自你上传的补充材料 ZIP：
- 41467_2023_37246_MOESM3_ESM (1).zip

## 这套包里包含什么
这是把“能用于建库的内容”整理成标准 CSV 的版本，重点保留：
1. 脑区字典与样本/脑区元数据
2. 表达基因目录
3. 各脑区 DEG 概览
4. WGCNA 模块注释
5. 模块功能、神经递质、细胞类型富集
6. CT 相关基因及其富集分析

## 这套包里不包含什么
这次补充材料里没有完整的 bulk atlas 主表达矩阵，因此这里仍然没有：
- gene × 819 samples
- gene × 110 regions

也就是说，这一包适合用于：
- 搭数据库骨架
- 建 region / gene / module / enrichment 注释层
- 为后续接入真正表达矩阵做准备

但还不能单独充当完整 bulk atlas reference matrix。

## 文件分组
- 01-03: 区域、脑区缩写、猴子信息
- 04-06: 表达基因目录与脑区 DEG 统计
- 07-18: 模块、通路、神经递质、细胞类型注释
- 19-31: CT 相关基因与富集
- 32-35: 我额外整理的建库辅助表

## 最推荐优先使用
- 02_region_dictionary.csv
- 33_region_annotation_joined.csv
- 32_region_sample_presence_matrix.csv
- 34_region_qc_summary.csv
- 35_recommended_database_tables.csv
