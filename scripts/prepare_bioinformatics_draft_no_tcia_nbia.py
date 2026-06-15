from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "manuscript" / "Bioinformatics_DraftA_bilingual_20260606.md"
OUTPUT = ROOT / "manuscript" / "Bioinformatics_DraftA_bilingual_no_TCIA_NBIA_20260609.md"


def replace_once(text: str, old: str, new: str) -> str:
    count = text.count(old)
    if count != 1:
        raise ValueError(f"Expected one occurrence, found {count}: {old[:80]!r}")
    return text.replace(old, new, 1)


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")

    replacements = [
        (
            "**Draft status:** Bilingual scientific-content draft; author details withheld; TCGA-TCIA MRI accuracy pending  ",
            "**Draft status:** Bilingual scientific-content draft; author details withheld  ",
        ),
        (
            "Glioma RNA-seq analyses further showed disease-associated prediction bias. Direct validation in clinical cfRNA cohorts and MRI-derived tumor-location validation remain future steps.",
            "Glioma RNA-seq analyses further showed disease-associated prediction bias. Direct validation in independent clinical cfRNA cohorts remains a future step.",
        ),
        (
            "胶质瘤RNA-seq进一步显示疾病相关预测偏倚。真实临床cfRNA队列及MRI肿瘤位置真值验证仍待完成。",
            "胶质瘤RNA-seq进一步显示疾病相关预测偏倚。真实独立临床cfRNA队列验证仍待完成。",
        ),
        (
            "Public datasets include the macaque atlas, AHBA, Ivy GAP, TCGA and TCIA, subject to their respective access conditions.",
            "Public datasets include the macaque atlas, AHBA, Ivy GAP and TCGA, subject to their respective access conditions.",
        ),
        (
            "所用公开资源包括猕猴图谱、AHBA、Ivy GAP、TCGA和TCIA，并遵循各数据库的访问条件。",
            "所用公开资源包括猕猴图谱、AHBA、Ivy GAP和TCGA，并遵循各数据库的访问条件。",
        ),
        (
            "### 2.8 Glioma RNA-seq and MRI linkage / 胶质瘤RNA-seq与MRI关联",
            "### 2.8 Glioma RNA-seq disease-domain analysis / 胶质瘤RNA-seq疾病域分析",
        ),
        (
            "\n**English.** TCGA patients were matched to TCIA MRI collections by patient identifier. Among 800 RNA-seq patients, 156 had an MRI collection match, 105 met minimum segmentation criteria and 73 had complete T1, contrast-enhanced T1, T2 and FLAIR imaging. Draft B will add DICOM conversion, multimodal quality control, tumor segmentation, atlas registration and MRI-derived Network labels before evaluating transcriptome-MRI agreement.\n\n**中文。** TCGA患者通过patient identifier与TCIA MRI collection匹配。800名RNA-seq患者中，156名匹配MRI collection，105名满足最低分割条件，73名具备完整T1、增强T1、T2和FLAIR。Draft B将在完成DICOM转换、多模态质控、肿瘤分割、atlas配准和MRI来源Network标签后，评价转录组与MRI的一致性。\n",
            "",
        ),
        (
            "Endpoint-specific evaluable sample counts were always reported. No TCGA location accuracy was calculated without MRI-derived ground truth.",
            "Endpoint-specific evaluable sample counts were always reported. Ivy GAP and TCGA disease-domain analyses summarized prediction distributions only and were not treated as anatomical accuracy analyses.",
        ),
        (
            "所有分析均报告终点特异的有效样本量。在没有MRI真值时不计算TCGA位置准确率。",
            "所有分析均报告终点特异的有效样本量。Ivy GAP和TCGA疾病域分析仅汇总预测分布，不作为解剖定位准确率分析。",
        ),
        (
            "Because Ivy labels represent tumor microenvironment structures and TCGA lacked MRI-derived location labels, these distributions cannot establish localization accuracy.",
            "Because Ivy labels represent tumor microenvironment structures and TCGA had no independent anatomical reference labels in this analysis, these distributions cannot establish localization accuracy.",
        ),
        (
            "由于Ivy标签代表肿瘤微环境结构，且TCGA尚缺MRI位置标签，这些分布不能用于证明定位准确率。",
            "由于Ivy标签代表肿瘤微环境结构，且本分析中的TCGA没有独立解剖参考标签，这些分布不能用于证明定位准确率。",
        ),
        (
            "\n### 3.7 Transcriptome-MRI validation cohort was established but not yet scored / 已建立转录组-MRI验证队列但尚未评分\n\n**English.** Of 800 TCGA RNA-seq patients, 156 matched TCIA MRI collections (19.5%). Seventy-three had complete four-modality MRI and 105 met minimum segmentation criteria. These cohorts define the independent validation path for Draft B. At the Draft A stage, all TCGA accuracy fields remain intentionally null because MRI-derived Network truth has not yet been generated.\n\n**中文。** 800名TCGA RNA-seq患者中，156名匹配TCIA MRI collection（19.5%）；其中73名具有完整四模态MRI，105名满足最低分割条件。这些队列构成Draft B的独立验证路径。在Draft A阶段，TCGA准确率字段有意保持为空，因为MRI来源的Network真值尚未生成。\n",
            "",
        ),
        (
            "Future analyses should mask tumor-, immune- and proliferation-associated genes, weight normal-brain-enriched genes, stratify GBM and LGG, and compare tumor center with edema-involved anatomy. MRI-derived labels are essential to determine whether any remaining signal corresponds to tumor location.",
            "Future analyses should mask tumor-, immune- and proliferation-associated genes, weight normal-brain-enriched genes, and stratify GBM and LGG. Independent anatomically labeled cohorts would be required to determine whether any remaining signal corresponds to tumor location.",
        ),
        (
            "后续应屏蔽肿瘤、免疫和增殖相关基因，提高正常脑富集基因权重，分层评价GBM与LGG，并区分肿瘤中心与水肿累及解剖。只有MRI来源标签才能判断剩余信号是否真正对应肿瘤位置。",
            "后续应屏蔽肿瘤、免疫和增殖相关基因，提高正常脑富集基因权重，并分层评价GBM与LGG。只有在独立且具有解剖标签的队列中，才能判断剩余信号是否真正对应肿瘤位置。",
        ),
        (
            "Third, MRI validation is incomplete, and no TCGA localization accuracy is reported. Fourth, the current workspace lacks commit-level public repository metadata. Finally, no direct cfRNA cohort was analyzed.",
            "Third, the unlabeled glioma analyses characterize domain shift but cannot establish localization accuracy. Fourth, the current workspace lacks commit-level public repository metadata. Finally, no direct cfRNA cohort was analyzed.",
        ),
        (
            "第三，MRI验证尚未完成，本文不报告TCGA定位准确率。第四，当前工作区缺少可公开引用的commit级仓库元数据。第五，尚未分析真实cfRNA队列。",
            "第三，无标签胶质瘤分析只能刻画域偏移，不能证明定位准确率。第四，当前工作区缺少可公开引用的commit级仓库元数据。第五，尚未分析真实cfRNA队列。",
        ),
        (
            "The framework provides a controlled foundation for future transcriptome-MRI validation and cfRNA source-tracing studies.",
            "The framework provides a controlled foundation for future cfRNA source-tracing studies and independent human validation.",
        ),
        (
            "该框架为后续转录组-MRI验证及cfRNA来源溯源研究提供了受控基础。",
            "该框架为后续cfRNA来源溯源研究及独立人源验证提供了受控基础。",
        ),
        (
            "**English.** Bo2023, AHBA, Ivy GAP, TCGA and TCIA data are available from their original repositories under the relevant access conditions. Exact accession identifiers, preprocessing manifests and restricted-access statements will be finalized before submission. TCIA restricted MRI data will not be redistributed.",
            "**English.** Bo2023, AHBA, Ivy GAP and TCGA data are available from their original repositories under the relevant access conditions. Exact accession identifiers, preprocessing manifests and access statements will be finalized before submission.",
        ),
        (
            "**中文。** Bo2023、AHBA、Ivy GAP、TCGA和TCIA数据可依据各原始数据库的访问条件获取。具体登录号、预处理清单和受限访问声明将在投稿前补全。受限TCIA MRI数据不会被再分发。",
            "**中文。** Bo2023、AHBA、Ivy GAP和TCGA数据可依据各原始数据库的访问条件获取。具体登录号、预处理清单和访问声明将在投稿前补全。",
        ),
        (
            "\n7. Clark K. et al. (2013) The Cancer Imaging Archive (TCIA): maintaining and operating a public information repository. *Journal of Digital Imaging*, 26, 1045-1057. doi:10.1007/s10278-013-9622-7.\n\n## Draft B placeholders / Draft B占位项\n\n1. MRI preprocessing and segmentation quality-control flowchart.\n2. MRI-derived Network truth definition.\n3. TCGA-TCIA Network Top1/Top3 accuracy.\n4. GBM versus LGG stratified performance.\n5. Tumor center versus whole-tumor/edema-involved region sensitivity analysis.\n6. Tumor/immune/proliferation gene masking ablation.\n7. Human-macaque ortholog and expression-scale calibration.\n",
            "\n\n## Future analyses / 后续分析\n\n1. GBM versus LGG stratified prediction-distribution analysis.\n2. Tumor-, immune- and proliferation-gene masking ablation.\n3. Human-macaque ortholog and expression-scale calibration.\n4. Independent anatomically labeled human validation.\n5. Direct validation in plasma or cerebrospinal-fluid cfRNA cohorts.\n",
        ),
    ]

    for old, new in replacements:
        text = replace_once(text, old, new)

    forbidden = ["TCIA", "NBIA", "MRI", "BraTS", "DICOM"]
    remaining = [term for term in forbidden if term.lower() in text.lower()]
    if remaining:
        raise ValueError(f"Forbidden terms remain in manuscript: {remaining}")

    OUTPUT.write_text(text, encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
