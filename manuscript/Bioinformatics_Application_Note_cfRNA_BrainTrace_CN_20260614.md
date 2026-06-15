# cfRNA-BrainTrace：基于灵长类转录组图谱的RNA-seq分层脑来源推断

**文章类型：** Application Note  
**栏目：** Gene expression  
**作者及单位：** [投稿前补充]  
**通讯作者：** [投稿前补充姓名及邮箱]

## 摘要

### 概要

cfRNA-BrainTrace是一套基于Python和Streamlit开发的应用程序，用于从组织或体液RNA表达谱中分层推断脑来源候选。软件首先将基因符号与灵长类脑参考图谱对齐，随后使用折内筛选的200基因Pearson相关模型，对10个宽粒度解剖Network进行评分；对于低margin样本，进一步执行pairwise rescue，并输出Network、broad anatomy、lobe及探索性的exact-region结果和置信度诊断。命令行和网页界面均可生成候选排序、审计表格及标准化结果摘要。在819个猕猴脑样本中，留一样本Network Top1和Top3准确率分别为55.8%和88.0%，留一猕猴准确率分别为53.2%和86.7%。正常人脑外部分析保留了部分粗粒度解剖信息；在配对胶质瘤转录组-MRI数据中，broad anatomy Top3候选覆盖率较高。对于缺少患者级解剖真值的体液RNA数据，软件明确将结果报告为迁移压力测试，而非定位准确率。cfRNA-BrainTrace因此提供了一套可重复实现，用于粗粒度脑来源候选排序及审计图谱驱动脑RNA溯源的分辨率边界。

### 可获得性与实现

运行环境为Python 3.11及以上版本，提供命令行和Streamlit界面。源代码、文档、测试、示例数据及带版本号的软件发布包将在以下地址公开：**[投稿前必须填写公开GitHub或永久归档地址]**。在线演示地址：**[投稿前必须填写稳定应用网址]**。软件许可证：**[投稿前必须填写OSI认可的开源许可证]**。

### 联系方式

**[投稿前必须填写通讯作者邮箱]**

### 补充材料

补充数据可在线获取。

## 1 引言

游离RNA中可保留可恢复的组织来源信息，但脑内解剖来源推断尤其困难，因为不同脑区共享转录程序，而体液RNA还具有稀释、稀疏及非脑组织混合等特征（Vorperian等，2022）。现有表达图谱方法通常直接返回单一标签，却不说明数据是否真正支持该解剖分辨率。在转化研究中，这会造成重要问题：宽粒度候选可能具有可重复性，但exact-region定位未必可靠。

cfRNA-BrainTrace实现了一种分层推断方案。软件使用Bo2023猕猴脑转录组图谱（Bo等，2023），从lobe、broad anatomy、10分类功能解剖Network及exact region四个层级对来源候选进行排序。应用程序将Network和粗粒度Top3候选作为主要输出，并同时显示低margin、低覆盖率及超出参考范围等情况。其目标用户是需要标准化评分接口、可重复诊断结果，以及在组织RNA-seq和前瞻性液体活检研究中明确结论边界的研究人员。

## 2 软件说明

### 2.1 输入与预处理

应用程序接受包含基因符号和TPM类丰度值的两列表达文件；可选的样本、受试者、诊断及解剖元数据会保留在导出报告中。软件对基因标识进行标准化，并与所选模型marker取交集。外部TPM类输入在评分前进行`log1p`转换。该处理属于跨尺度相关分析，并非DESeq2 variance-stabilized数值的重建。软件同时报告匹配marker数量及非零覆盖率，从而区分稀疏体液输入与marker覆盖充分的组织表达谱。

### 2.2 分层推断

锁定模型基于819个猕猴脑样本建立，覆盖9只猕猴和110个脑区。在每个验证折内独立筛选200个判别基因，随后使用Pearson相关计算样本与各类别centroid之间的相似性。当Top1与Top2相关分数的margin不高于0.002时，pairwise rescue模型会重新评价排名前三的Network。Rescue阈值与模型共同保存，并由命令行和网页界面以相同方式调用。

预测结果通过带版本号的解剖字典进一步映射到lobe及broad anatomy。Exact-region结果属于探索性输出，并受到候选Network约束，不作为独立的临床终点。当前10分类Network参考不包含小脑，因此小脑或后颅窝样本会被标记为超出参考范围，而不会强制映射至最接近的现有类别。

### 2.3 界面与输出

Streamlit界面支持文件上传、元数据检查、Network及分层候选展示、结果表下载和模型警告。命令行界面支持自动化评分、benchmark导出及可重复结果打包。每个样本报告包含Top1和Top3候选、相关分数、经softmax归一化的展示概率、Top1-Top2 margin、标准化熵、marker覆盖率、rescue状态和各解剖层级标签。批量导出文件同时保存模型版本、输入设置及作图源数据表。

应用程序将预测与验证严格区分。只有提供独立定义的解剖真值标签时，软件才计算准确率。对于无真值样本，仅报告候选分布、置信度、熵、稳定性和marker覆盖率，避免将无标签体液队列描述为定位验证。

## 3 软件实现

cfRNA-BrainTrace使用Python开发，采用NumPy、pandas和SciPy完成评分，使用scikit-learn相关工具完成评价，使用Plotly和Matplotlib进行可视化，并以Streamlit构建网页界面。模型以带版本号的NumPy和JSON文件发布，其中包含marker顺序、类别centroid、解剖映射和rescue参数。正式Network路线实现在`core/network_tracing.py`中；命令行和网页应用调用相同的核心评分函数，以避免不同界面产生不一致结果。

代码仓库包含模型加载、marker对齐、非有限值处理、概率归一化及确定性排序等单元测试。Benchmark脚本可复现源域验证，并同时导出图件及其CSV源数据。探索性域适配代码被单独保留，默认评分器不会调用该路线，因为其性能依赖评价终点，并且部分步骤使用了目标队列信息，详见补充方法。

## 4 验证结果

严格留一样本验证的Network Top1准确率为55.8%，Top3准确率为88.0%。留一猕猴验证的Top1和Top3准确率分别为53.2%和86.7%，说明主要候选排序性能并非仅由同一个体内重复取样造成（图1B）。

跨物种评价采用Allen Human Brain Atlas的协调标签（Hawrylycz等，2012）。在233个具有支持性粗粒度标签的样本中，Network Top1和Top3准确率分别为32.6%和55.4%，lobe Top1准确率为44.2%。更细粒度的跨物种迁移性能明显较弱，支持软件采用分层结果报告策略。

在65例具有配对TCGA-LGG RNA-seq和BraTS肿瘤分割的患者中（Bakas等，2017），broad anatomy Top3 strict和tolerant覆盖率分别为75.4%和83.1%；在64例参考范围内患者中，Network Top3分别为21.9%和35.9%（图1C）。这些结果支持在肿瘤组织中开展粗粒度候选排序，但不支持单一脑区定位。

三套公开液体活检队列分别包含血清细胞外囊泡RNA、肿瘤相关细胞外囊泡RNA及脑脊液RNA，均被作为外部迁移压力测试。尽管独立的TPM与counts重算CPM审计中67项数值实现检查全部通过，这些队列的预测通常仍表现为低margin、高熵或预处理敏感。由于缺少患者级影像真值，本研究未计算其定位准确率。完整验证设计、置信区间、疾病域分析及adapted路线的负结果见补充材料。

## 5 结论

cfRNA-BrainTrace将经过验证的图谱相关分析流程封装为可复用的命令行和网页应用，并将解剖分辨率作为预测结果的明确组成部分。其主要用途是提供粗粒度脑来源候选排序，以及透明的置信度和marker覆盖诊断。该软件不会把无标签体液RNA的预测转化为临床定位结论。未来版本仍需建立人源体液背景模型、校准拒绝规则，并在具有患者级解剖真值的独立队列中完成验证。

## 基金资助

**[投稿前补充基金信息]**

## 利益冲突

作者声明不存在利益冲突。

## 数据可获得性

Bo2023猕猴转录组图谱、Allen Human Brain Atlas、TCGA、BraTS、Ivy Glioblastoma Atlas Project及GEO数据均可根据原始数据资源规定获取，具体登录号和访问条件见补充材料。复现本文结果所需的处理后评价表将包含在带版本号的软件发布包中。

## 参考文献

Bakas,S.等（2017）Advancing The Cancer Genome Atlas glioma MRI collections with expert segmentation labels and radiomic features. *Scientific Data*, **4**, 170117.

Bo,T.等（2023）Brain-wide and cell-specific transcriptomic insights into MRI-derived cortical morphology in macaque monkeys. *Nature Communications*, **14**, 1283.

Hawrylycz,M.J.等（2012）An anatomically comprehensive atlas of the adult human brain transcriptome. *Nature*, **489**, 391-399.

Vorperian,S.K.等（2022）Cell types of origin of the cell-free transcriptome. *Nature Biotechnology*, **40**, 855-861.

## 图注

**图1. cfRNA-BrainTrace工作流程与验证。**（A）软件输入、锁定评分路线及分层输出。命令行和Streamlit界面调用相同的核心评分模块。（B）严格留一样本（LOSO）及留一猕猴（LOMO）验证中的Network Top1和Top3准确率。（C）正常人脑及配对TCGA-LGG/BraTS数据中的外部粗粒度性能。AHBA报告协调后的Network准确率；配对胶质瘤报告strict及tolerant Top3覆盖率。（D）结果解释原则：具有独立真值的组织数据可计算准确率；无标签EV-RNA及CSF-RNA仅用于迁移、置信度和稳定性诊断。

## 表格

**表1. 软件主要功能与输出。**

| 组成部分 | 实现方式 | 主要输出 |
| --- | --- | --- |
| 输入 | 基因符号表达表及可选元数据 | Marker匹配数量及非零覆盖率 |
| 评分 | 200基因Pearson相关加pairwise rescue | Network Top1/Top3排序及rescue状态 |
| 分层结构 | Lobe、broad anatomy、Network及exact region | 各分辨率候选及超出参考范围标志 |
| 诊断 | 相关分数、展示概率、margin及熵 | 置信度及迁移质量审计 |
| 界面 | Streamlit网页应用和命令行界面 | 交互式检查及可重复批量导出 |
| 验证原则 | 仅在具有独立解剖真值时计算准确率 | 无标签体液样本报告为压力测试 |
