# Bo2023 buildkit 已整合进 SQLite

已将 `bo2023_bulk_atlas_buildkit` 全部 CSV 导入 `cfrna_source_tracing.db`。

## 新增 SQLite 元数据表
- bo2023_buildkit_catalog
- bo2023_buildkit_column_map

## 新增关键视图
- bo2023_regions
- bo2023_region_qc_overview
- bo2023_region_sample_coverage
- bo2023_ct_gene_catalog

## 可直接查询示例
```sql
SELECT * FROM bo2023_regions LIMIT 10;
SELECT * FROM bo2023_region_qc_overview ORDER BY rin DESC LIMIT 10;
SELECT region, total_sample_count FROM bo2023_region_sample_coverage ORDER BY total_sample_count DESC;
SELECT * FROM bo2023_buildkit_catalog ORDER BY table_name;
```

## 说明
这些数据是 bulk atlas 的建库骨架与注释层，尚不包含完整 `gene × region` 表达矩阵，因此不会自动替换当前 legacy tracing reference。


## 2026-04-01 前端修复
- 修复 Atlas 浏览器脑区总览查询错误：`full_name` 等字段实际位于 `bo2023_regions` 视图，而不是 `bo2023_region_sample_coverage`。
- 现已改为从 `bo2023_regions`、`bo2023_region_sample_coverage`、`bo2023_region_qc_overview` 三个视图联合查询。
