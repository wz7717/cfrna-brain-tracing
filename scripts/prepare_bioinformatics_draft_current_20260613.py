from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "manuscript" / "Bioinformatics_DraftA_bilingual_no_TCIA_NBIA_20260609.md"
OUTPUT = ROOT / "manuscript" / "Bioinformatics_DraftB_candidate_bilingual_20260613.md"


def replace_once(text: str, old: str, new: str) -> str:
    count = text.count(old)
    if count != 1:
        raise ValueError(f"Expected one occurrence, found {count}: {old[:100]!r}")
    return text.replace(old, new, 1)


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    replacements = [
        (
            "**Draft status:** Bilingual scientific-content draft; author details withheld  ",
            "**Draft status:** Bilingual Draft B candidate; author details withheld; updated through 13 June 2026  ",
        ),
        (
            "**English.** The Network model used 200 fold-selected discriminative genes and Pearson correlation, followed by pairwise correlation rescue for low-margin Top1 predictions. In strict leave-one-sample-out validation of 819 macaque brain samples, Network Top1 accuracy increased from 52.7% to 55.8%, while Top3 accuracy remained 88.0% (45 gains, 20 losses, paired P=0.0026). Leave-one-monkey-out validation yielded 53.2% Top1 and 86.7% Top3 accuracy. Region-group and exact-region performance was lower and was treated as secondary and exploratory, respectively. In label-harmonized Allen Human Brain Atlas RNA-seq samples, coarse Network Top1/Top3 accuracy was 24.9%/55.4%, demonstrating transferable signal but substantial cross-species domain shift. Glioma RNA-seq analyses further showed disease-associated prediction bias. Direct validation in independent clinical cfRNA cohorts remains a future step.",
            "**English.** The Network model used 200 fold-selected discriminative genes and Pearson correlation, followed by pairwise correlation rescue for low-margin Top1 predictions. In strict leave-one-sample-out validation of 819 macaque brain samples, Network Top1 accuracy increased from 52.7% to 55.8%, while Top3 accuracy remained 88.0% (paired P=0.0026). Leave-one-monkey-out validation yielded 53.2% Top1 and 86.7% Top3 accuracy. In 233 label-harmonized Allen Human Brain Atlas samples, coarse Network Top1/Top3 accuracy was 32.6%/55.4% and lobe Top1 was 44.2%. In 65 patients with paired TCGA-LGG RNA-seq and BraTS tumor segmentations, the locked baseline achieved lobe Top3 strict/tolerant coverage of 84.6%/89.2% and broad-anatomy Top3 of 75.4%/83.1%, whereas Network Top3 was 21.9%/35.9%. These results support coarse candidate ranking rather than reliable single-point tumor localization. Direct validation in independent clinical cfRNA cohorts remains a future step.",
        ),
        (
            "**中文。** Network模型采用折内筛选的200个判别基因及Pearson相关性，并对Top1分差较小的样本启动pairwise correlation rescue。在819个猕猴脑样本的严格留一样本验证中，Network Top1准确率由52.7%提高至55.8%，Top3保持88.0%（45个增益、20个损失，配对检验P=0.0026）。留一猴验证的Top1和Top3分别为53.2%和86.7%。Region Group与Exact Region性能较低，分别作为次级和探索性结果。经标签协调的Allen Human Brain Atlas人脑RNA-seq中，Network coarse Top1/Top3为24.9%/55.4%，提示存在可迁移信号，同时也显示明显跨物种域偏移。胶质瘤RNA-seq进一步显示疾病相关预测偏倚。真实独立临床cfRNA队列验证仍待完成。",
            "**中文。** Network模型采用折内筛选的200个判别基因及Pearson相关性，并对Top1分差较小的样本启动pairwise correlation rescue。在819个猕猴脑样本的严格留一样本验证中，Network Top1准确率由52.7%提高至55.8%，Top3保持88.0%（配对P=0.0026）；留一猴验证的Top1和Top3分别为53.2%和86.7%。在233个标签协调的Allen Human Brain Atlas样本中，coarse Network Top1/Top3为32.6%/55.4%，lobe Top1为44.2%。在65例具有配对TCGA-LGG RNA-seq和BraTS肿瘤分割的患者中，锁定baseline的lobe Top3 strict/tolerant为84.6%/89.2%，broad anatomy Top3为75.4%/83.1%，而Network Top3仅为21.9%/35.9%。结果支持粗粒度候选区域排序，不支持可靠的单点肿瘤定位。真实独立临床cfRNA队列验证仍待完成。",
        ),
        (
            "Public datasets include the macaque atlas, AHBA, Ivy GAP and TCGA, subject to their respective access conditions.",
            "Public datasets include the macaque atlas, AHBA, Ivy GAP, TCGA and the paired BraTS glioma imaging resource, subject to their respective access conditions.",
        ),
        (
            "所用公开资源包括猕猴图谱、AHBA、Ivy GAP和TCGA，并遵循各数据库的访问条件。",
            "所用公开资源包括猕猴图谱、AHBA、Ivy GAP、TCGA及配对的BraTS胶质瘤影像资源，并遵循各数据库的访问条件。",
        ),
        (
            "We then examined transfer to AHBA normal human brain and characterized domain shift in Ivy GAP and TCGA glioma RNA-seq.",
            "We then examined transfer to AHBA normal human brain, evaluated a locked baseline in patients with paired TCGA-LGG RNA-seq and BraTS tumor segmentations, and characterized domain shift in Ivy GAP and the full TCGA glioma cohorts.",
        ),
        (
            "随后在AHBA正常人脑中评价跨物种迁移，并在Ivy GAP和TCGA胶质瘤RNA-seq中刻画域偏移。",
            "随后在AHBA正常人脑中评价跨物种迁移，在配对TCGA-LGG RNA-seq与BraTS肿瘤分割患者中评价锁定baseline，并在Ivy GAP及完整TCGA胶质瘤队列中刻画域偏移。",
        ),
        (
            "### 2.8 Glioma RNA-seq disease-domain analysis / 胶质瘤RNA-seq疾病域分析\n\n**English.** Ivy GAP comprised 122 samples from five glioblastoma microanatomical structures. These labels describe tumor compartments rather than normal anatomical origin; therefore, only prediction distributions were summarized. TCGA-GBM/LGG expression processing yielded 801 RNA-seq samples from 800 patients. Predictions were generated, but accuracy was not computed because anatomical truth was unavailable.\n\n**中文。** Ivy GAP包含122个样本和5类胶质母细胞瘤微解剖结构。这些标签描述肿瘤区室而非正常脑区来源，因此仅汇总预测分布。TCGA-GBM/LGG表达数据共得到801个RNA-seq样本、800名患者。模型已生成预测，但由于缺少解剖真值，未计算准确率。\n\n### 2.9 Statistical analysis / 统计分析\n\n**English.** The primary metrics were sample-level Top1 and Top3 accuracy for Network. Secondary metrics were corresponding Region Group values, and Exact Region metrics were exploratory. Paired gains and losses were calculated by comparing baseline and rescued predictions on the same held-out samples. Endpoint-specific evaluable sample counts were always reported. Ivy GAP and TCGA disease-domain analyses summarized prediction distributions only and were not treated as anatomical accuracy analyses.\n\n**中文。** 主要指标为样本级Network Top1和Top3准确率；Region Group对应指标为次级结果，Exact Region为探索性结果。配对增益与损失通过同一留出样本的基线和rescue预测比较计算。所有分析均报告终点特异的有效样本量。Ivy GAP和TCGA疾病域分析仅汇总预测分布，不作为解剖定位准确率分析。",
            """### 2.8 Glioma RNA-seq disease-domain analysis / 胶质瘤RNA-seq疾病域分析

**English.** Ivy GAP comprised 122 samples from five glioblastoma microanatomical structures. These labels describe tumor compartments rather than normal anatomical origin; therefore, only prediction distributions were summarized. The full TCGA cohort contained 285 GBM samples from 284 patients and 516 LGG samples from 516 patients. GBM-LGG differences in patient-level Top1 Network distributions were assessed by permutation testing, Cramer's V and Jensen-Shannon divergence. Low confidence was prespecified as a Top1-Top2 score margin below 0.002.

**中文。** Ivy GAP包含122个样本和5类胶质母细胞瘤微解剖结构；这些标签描述肿瘤区室而非正常脑区来源，因此仅汇总预测分布。完整TCGA队列包含285个GBM样本（284名患者）和516个LGG样本（516名患者）。患者级Top1 Network分布的GBM-LGG差异采用置换检验、Cramer's V和Jensen-Shannon divergence评价；Top1-Top2得分margin低于0.002预设为低置信度。

### 2.9 Paired TCGA-LGG transcriptome and BraTS validation / TCGA-LGG转录组与BraTS配对验证

**English.** Sixty-five TCGA-LGG patients had paired RNA-seq, pre-operative multimodal MRI and tumor segmentations; 62 segmentations were manually corrected and three were automatic. Tumor masks and the TZO116 atlas were checked in the same image space. Anatomical truth used direct tumor-mask overlap with valid atlas labels without unbounded nearest-label filling. Median direct atlas coverage was 66.7%; two patients had coverage below 50%. One cerebellar tumor was outside the current 10-Network reference, leaving 64 Network-evaluable patients.

**中文。** 65名TCGA-LGG患者具有配对RNA-seq、术前多模态MRI及肿瘤分割，其中62例为人工校正分割、3例为自动分割。肿瘤mask与TZO116 atlas在同一影像空间内完成核查。解剖真值仅使用肿瘤mask与有效atlas标签的直接重叠，不采用无界最近邻补标签。直接atlas覆盖率中位数为66.7%，2例低于50%；1例小脑肿瘤超出当前10类Network参考范围，因此Network有效样本为64例。

**English.** Strict accuracy required agreement with the dominant directly overlapping label. Tolerant accuracy accepted any directly overlapping candidate label meeting the prespecified overlap rule. The locked baseline was evaluated at lobe, broad-anatomy and Network levels using Top1 and Top3. Exact-region accuracy was not evaluated in this cohort.

**中文。** Strict准确率要求预测与直接重叠的主导标签一致；tolerant准确率允许命中符合预设重叠规则的任一直接重叠候选标签。锁定baseline在lobe、broad anatomy和Network层级评价Top1及Top3；本队列不评价Exact Region准确率。

### 2.10 Domain-adaptation sensitivity analysis / 域适配敏感性分析

**English.** Tumor-adapted routes used the target cohort's unlabeled expression distribution for gene filtering and scale harmonization. Raw-TPM, log1p and harmonized variants were compared with the locked baseline. A calibrated harmonized route used the 65 MRI labels to select calibration parameters within nested leave-one-out analysis and was therefore classified as cohort-internal calibration rather than independent external validation. A source-derived out-of-distribution threshold was also tested.

**中文。** 肿瘤适配路线使用目标队列的无标签表达分布进行基因过滤和尺度协调，并比较raw TPM、log1p及harmonized变体与锁定baseline。calibrated harmonized路线在嵌套留一分析中使用65例MRI标签选择校准参数，因此属于队列内校准，而非独立外部验证。此外还测试了来源于猕猴数据的out-of-distribution阈值。

### 2.11 Statistical analysis / 统计分析

**English.** The primary internal metrics were sample-level Network Top1 and Top3 accuracy. Region Group was secondary and Exact Region exploratory. In the paired TCGA-LGG/BraTS cohort, the principal interpretable disease-domain readout was broad-anatomy Top3 strict/tolerant coverage, with lobe Top3 supportive and Network metrics secondary. Binomial proportions were reported with 95% confidence intervals. Ivy GAP and unpaired TCGA analyses were treated as stress tests rather than localization-accuracy analyses.

**中文。** 内部验证的主要指标为样本级Network Top1和Top3；Region Group为次级终点，Exact Region为探索性终点。在配对TCGA-LGG/BraTS队列中，疾病域主要可解释读出为broad anatomy Top3 strict/tolerant覆盖率，lobe Top3为支持性结果，Network指标为次要结果。二项比例报告95%置信区间。Ivy GAP及无配对真值的TCGA分析被视为压力测试，而非定位准确率分析。""",
        ),
        (
            "**English.** Among 233 AHBA samples with harmonized coarse labels, Network Top1 accuracy was 24.9% and Top3 accuracy was 55.4%. Broad-anatomy/lobe Top1 accuracy was 44.2%. Among 91 exact-mapped samples, Exact Region Top1 and Top3 were 9.9% and 29.7%. The results show that a macaque-derived model retains coarse anatomical information in human normal-brain RNA-seq but undergoes a marked reduction relative to within-species validation.",
            "**English.** Among 233 AHBA samples with harmonized coarse labels, Network Top1 accuracy was 32.6% and Top3 accuracy was 55.4%. Lobe Top1 accuracy was 44.2%. Among 91 exact-mapped samples, Exact Region Top1 and Top3 were 9.9% and 29.7%. Moreover, 92.2% of Top1 outputs were designated low-resolution or recommended for manual review. The results show that a macaque-derived model retains coarse anatomical information in normal human brain RNA-seq, while fine localization remains unstable.",
        ),
        (
            "**中文。** 在233个具有协调coarse标签的AHBA样本中，Network Top1为24.9%，Top3为55.4%；宽解剖/脑叶Top1为44.2%。在91个exact-mapped样本中，Exact Region Top1/Top3为9.9%/29.7%。结果表明，猕猴来源模型在人正常脑RNA-seq中仍保留宽粒度解剖信息，但相较同物种内部验证明显下降。",
            "**中文。** 在233个具有协调coarse标签的AHBA样本中，Network Top1为32.6%，Top3为55.4%，lobe Top1为44.2%。在91个exact-mapped样本中，Exact Region Top1/Top3为9.9%/29.7%。此外，92.2%的Top1输出被标记为低分辨率或建议人工复核。结果表明，猕猴来源模型在人正常脑RNA-seq中仍保留宽粒度解剖信息，但精细定位仍不稳定。",
        ),
        (
            "### 3.6 Glioma expression exposed disease-domain shift / 胶质瘤表达揭示疾病域偏移\n\n**English.** Ivy GAP and TCGA glioma predictions were enriched for a limited subset of deep-brain and hippocampal-related outputs. Because Ivy labels represent tumor microenvironment structures and TCGA had no independent anatomical reference labels in this analysis, these distributions cannot establish localization accuracy. Instead, they indicate that tumor, immune, proliferative and cellular-composition signals may dominate similarity to a normal macaque reference.\n\n**中文。** Ivy GAP和TCGA胶质瘤预测集中于少数深部脑区及海马相关输出。由于Ivy标签代表肿瘤微环境结构，且本分析中的TCGA没有独立解剖参考标签，这些分布不能用于证明定位准确率。相反，它们提示肿瘤、免疫、增殖及细胞组成信号可能主导其与正常猕猴参考图谱的相似性。",
            """### 3.6 GBM and LGG showed measurable disease-domain differences / GBM与LGG存在可测量的疾病域差异

**English.** Patient-level Top1 Network distributions differed between 284 GBM and 516 LGG patients (permutation P=0.000600; Cramer's V=0.158; Jensen-Shannon divergence=0.0266 bits). Low-confidence predictions occurred in 38.4% of GBM and 26.0% of LGG patients (risk ratio 1.776). The effect size was modest, but GBM had consistently smaller probability and raw-score margins. Ivy GAP and unpaired TCGA outputs remain stress-test distributions and do not establish localization accuracy.

**中文。** 284名GBM与516名LGG患者的Top1 Network分布存在差异（置换P=0.000600；Cramer's V=0.158；Jensen-Shannon divergence=0.0266 bits）。低置信度预测在GBM和LGG中分别占38.4%和26.0%（风险比1.776）。效应量较小，但GBM的概率margin和原始得分margin持续更低。Ivy GAP及无配对真值的TCGA输出仍属于压力测试分布，不能证明定位准确率。

### 3.7 Paired tumor transcriptome-MRI validation favored coarse Top3 coverage / 配对肿瘤转录组-MRI验证支持粗粒度Top3覆盖

**English.** In 65 paired TCGA-LGG/BraTS patients, the locked baseline achieved lobe Top1 strict/tolerant accuracy of 9.2%/12.3% and Top3 of 84.6%/89.2%. Broad-anatomy Top1 was 6.2%/10.8% and Top3 was 75.4%/83.1%. Among 64 Network-evaluable patients, Network Top1 was 3.1%/4.7% and Top3 was 21.9%/35.9%. Thus, candidate coverage was strong at lobe and broad-anatomy Top3 but did not translate into reliable Top1 or fine Network localization.

**中文。** 在65例配对TCGA-LGG/BraTS患者中，锁定baseline的lobe Top1 strict/tolerant为9.2%/12.3%，Top3为84.6%/89.2%；broad anatomy Top1为6.2%/10.8%，Top3为75.4%/83.1%。在64例Network可评价患者中，Network Top1为3.1%/4.7%，Top3为21.9%/35.9%。因此，lobe和broad anatomy层级具有较强的Top3候选覆盖，但不能转化为可靠的Top1或精细Network定位。

### 3.8 Domain adaptation redistributed the Top1-Top3 trade-off / 域适配改变Top1与Top3权衡

**English.** Unsupervised harmonization increased tolerant Network Top1 from 4.7% to 20.3% and tolerant Network Top3 from 35.9% to 48.4%, but reduced tolerant broad-anatomy Top3 from 83.1% to 60.0%. Cohort-calibrated harmonization yielded 20.3% Network Top1, 46.9% Network Top3 and 70.8% broad-anatomy Top3. Because the first route used the target expression distribution and the second additionally used cohort labels, these analyses are sensitivity analyses rather than locked independent validation. The macaque-derived OOD threshold rejected all 65 patients.

**中文。** 无监督harmonization将tolerant Network Top1由4.7%提高至20.3%，Network Top3由35.9%提高至48.4%，但将tolerant broad anatomy Top3由83.1%降低至60.0%。队列内校准的harmonized路线获得20.3%的Network Top1、46.9%的Network Top3和70.8%的broad anatomy Top3。由于前者使用了目标表达分布，后者还使用了队列标签，这些结果属于敏感性分析而非锁定的独立外部验证。来源于猕猴的OOD阈值拒绝了全部65例患者。""",
        ),
        (
            "**English.** The glioma analyses emphasize that anatomical source and disease state are not interchangeable. Bulk tumor RNA contains malignant, immune, vascular and stromal components, and these programs may overwhelm residual regional identity. Future analyses should mask tumor-, immune- and proliferation-associated genes, weight normal-brain-enriched genes, and stratify GBM and LGG. Independent anatomically labeled cohorts would be required to determine whether any remaining signal corresponds to tumor location.",
            "**English.** The paired TCGA-LGG/BraTS results refine the disease-domain interpretation. High lobe and broad-anatomy Top3 coverage indicates that the locked model can retain coarse candidate information in bulk tumor RNA, but very low Top1 and Network performance preclude single-site localization claims. The GBM-LGG comparison and complete rejection by the macaque-derived OOD threshold further show that disease state and species shift materially alter confidence. Unsupervised harmonization improved selected Network metrics at the cost of broad Top3 coverage, illustrating that domain adaptation redistributes rather than uniformly resolves error.",
        ),
        (
            "**中文。** 胶质瘤分析强调了解剖来源与疾病状态不可互换。Bulk肿瘤RNA包含恶性、免疫、血管及基质成分，这些程序可能掩盖残余区域身份。后续应屏蔽肿瘤、免疫和增殖相关基因，提高正常脑富集基因权重，并分层评价GBM与LGG。只有在独立且具有解剖标签的队列中，才能判断剩余信号是否真正对应肿瘤位置。",
            "**中文。** 配对TCGA-LGG/BraTS结果进一步限定了疾病域解释。较高的lobe和broad anatomy Top3覆盖表明，锁定模型可在bulk肿瘤RNA中保留粗粒度候选信息，但极低的Top1和Network性能不支持单点定位。GBM-LGG差异及猕猴来源OOD阈值拒绝全部患者进一步说明疾病状态和物种偏移会显著改变模型置信度。无监督harmonization改善了部分Network指标，却损失broad Top3覆盖，说明域适配重新分配误差，而非全面解决误差。",
        ),
        (
            "Third, the unlabeled glioma analyses characterize domain shift but cannot establish localization accuracy.",
            "Third, the paired glioma cohort was limited to 65 LGG patients, two had low direct atlas coverage, and one cerebellar tumor was outside the reference Network space; strict and tolerant rules therefore require joint reporting. Fourth, transductive and cohort-calibrated adaptation cannot be interpreted as locked independent validation.",
        ),
        (
            "第三，无标签胶质瘤分析只能刻画域偏移，不能证明定位准确率。第四，当前工作区缺少可公开引用的commit级仓库元数据。第五，尚未分析真实cfRNA队列。",
            "第三，配对胶质瘤队列仅包含65例LGG患者，其中2例直接atlas覆盖率较低，1例小脑肿瘤超出参考Network范围，因此必须同时报告strict与tolerant规则。第四，transductive及队列内校准的适配结果不能解释为锁定独立验证。第五，当前工作区缺少可公开引用的commit级仓库元数据。第六，尚未分析真实cfRNA队列。",
        ),
        (
            "The framework provides a controlled foundation for future cfRNA source-tracing studies and independent human validation.",
            "The framework supports coarse candidate ranking across normal-human and paired glioma settings, while defining the calibration and validation required before cfRNA source tracing.",
        ),
        (
            "该框架为后续cfRNA来源溯源研究及独立人源验证提供了受控基础。",
            "该框架支持在正常人脑及配对胶质瘤场景中进行粗粒度候选排序，同时明确了进入cfRNA来源溯源前仍需完成的校准与验证。",
        ),
        (
            "Bo2023, AHBA, Ivy GAP and TCGA data are available from their original repositories under the relevant access conditions.",
            "Bo2023, AHBA, Ivy GAP, TCGA and BraTS data are available from their original repositories under the relevant access conditions.",
        ),
        (
            "Bo2023、AHBA、Ivy GAP和TCGA数据可依据各原始数据库的访问条件获取。",
            "Bo2023、AHBA、Ivy GAP、TCGA及BraTS数据可依据各原始数据库的访问条件获取。",
        ),
        (
            "## Future analyses / 后续分析\n\n1. GBM versus LGG stratified prediction-distribution analysis.\n2. Tumor-, immune- and proliferation-gene masking ablation.\n3. Human-macaque ortholog and expression-scale calibration.\n4. Independent anatomically labeled human validation.\n5. Direct validation in plasma or cerebrospinal-fluid cfRNA cohorts.",
            "## Future analyses / 后续分析\n\n1. Predefined tumor-, immune- and proliferation-gene masking ablation with source-domain protection metrics.\n2. Human-macaque ortholog and expression-scale calibration in donor-isolated validation.\n3. External replication of the paired transcriptome-anatomy analysis in an independent glioma cohort.\n4. Prospective validation of coarse candidate ranking in independently labeled human samples.\n5. Direct validation in plasma or cerebrospinal-fluid cfRNA cohorts.",
        ),
    ]
    for old, new in replacements:
        text = replace_once(text, old, new)
    text = text.replace(
        "Fourth, the current workspace lacks commit-level public repository metadata. Finally, no direct cfRNA cohort was analyzed.",
        "Fifth, the current workspace lacks commit-level public repository metadata. Sixth, no direct cfRNA cohort was analyzed.",
    )
    OUTPUT.write_text(text, encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
