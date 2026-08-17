from __future__ import annotations

import streamlit as st

from app.i18n import tr
from app.shared import render_page_hero
from core.model_lock import EXPECTED_MODEL_LOCK_ID


def display_help_page() -> None:
    render_page_hero(
        tr("帮助与教程", "Help & Tutorial"),
        tr(
            "在五分钟内运行 BrainTrace，并正确解释三层候选结果。",
            "Run BrainTrace in about five minutes and interpret its three candidate tiers correctly.",
        ),
        eyebrow="BrainTrace v0.1.16",
        pills=["Network Top3", "Resolution Group Top3", "Exact-Region exploratory Top3"],
    )
    st.header(tr("BrainTrace 是什么？", "What is BrainTrace?"))
    st.write(
        tr(
            "BrainTrace 使用锁定的灵长类转录组参考，将 RNA 表达谱排序为分层的脑来源候选。",
            "BrainTrace ranks RNA expression profiles into hierarchical brain-origin candidates using a locked primate transcriptomic reference.",
        )
    )

    st.header(tr("预期用途", "Intended use"))
    st.write(
        tr(
            "用于研究性候选排序和方法学探索；不是诊断工具，也不产生确定性解剖定位结论。",
            "Research candidate ranking and methodological exploration only; it is not diagnostic and does not produce deterministic anatomical localization calls.",
        )
    )

    st.header(tr("快速开始", "Quick Start"))
    st.markdown(
        tr(
            "1. 上传 raw-count 或 logCPM 表格。\n2. 检查模型基因覆盖率。\n3. 运行 BrainTrace。\n4. 首先阅读 Network Top3。\n5. 将 exact-region 输出仅作为探索性结果。",
            "1. Upload a raw-count or logCPM table.\n2. Check model-gene coverage.\n3. Run BrainTrace.\n4. Read Network Top3 first.\n5. Treat exact-region output as exploratory.",
        )
    )

    st.header(tr("输入格式", "Input format"))
    st.code("gene_symbol\traw_counts\nSATB2\t100", language="text")
    st.code("gene_symbol\tlogCPM\nSATB2\t4.25", language="text")
    st.info(
        tr(
            "推荐 raw counts 或 logCPM。TPM/logTPM 仅作为 fallback/兼容输入，不与正式路线等价。",
            "Raw counts or logCPM are recommended. TPM/logTPM is accepted only as fallback/compatibility input and is not equivalent to the formal route.",
        )
    )

    st.header(tr("加载示例数据", "Load example data"))
    st.write(
        tr(
            "在溯源分析页点击“加载示例数据”，再点击“运行示例”。示例不会进入样本管理或写入 SQLite。",
            "On the Tracing page, click Load example data and then Run example. The example does not enter Sample Management or write to SQLite.",
        )
    )
    st.warning(
        tr(
            "合成软件示例，仅用于演示软件运行，不代表生物学验证数据。",
            "Synthetic software example. Not biological validation data.",
        )
    )

    st.header(tr("三层结果解释", "Interpretation of the three tiers"))
    st.markdown(
        tr(
            "- **Network Top3**：主要、已验证的候选层级。\n- **Resolution Group Top3**：推荐报告的可分辨候选层级。\n- **Exact Region Top3**：仅探索性结果；不是确定性定位。",
            "- **Network Top3**: primary validated candidate tier.\n- **Resolution Group Top3**: recommended resolvable candidate tier.\n- **Exact Region Top3**: exploratory only; it is not a deterministic localization call.",
        )
    )

    st.header(tr("输入覆盖率要求", "Input coverage requirements"))
    st.markdown(
        tr(
            "Network 推理至少需要锁定 200 基因面板的 50%（≥100/200）。精细层级至少需要与区域参考重叠 20 个基因。低于门槛时返回 insufficient traceability，不强制预测。",
            "Network inference requires at least 50% of the locked 200-gene panel (≥100/200). Fine-tier tracing requires at least 20 genes overlapping the regional reference. Below either gate BrainTrace returns insufficient traceability rather than forcing a prediction.",
        )
    )

    st.header(tr("置信度诊断", "Confidence diagnostics"))
    st.write(
        tr(
            "结合 coverage、entropy 和 Top1–Top2 score margin 判断输入支持度与候选分离程度；这些指标不把候选排名转化为确定性定位。",
            "Use coverage, entropy and the Top1–Top2 score margin to assess input support and candidate separation; these diagnostics do not turn a ranking into deterministic localization.",
        )
    )

    st.header(tr("局限性", "Limitations"))
    st.warning(
        tr(
            "缺少患者级解剖真值的生物流体数据不能解释为定位准确率验证。跨物种、疾病和测序平台迁移均可能产生 domain shift。",
            "Biofluid datasets without patient-level anatomical truth must not be interpreted as localization-accuracy validation. Cross-species, disease and assay transfer may introduce domain shift.",
        )
    )

    st.header(tr("可复现性", "Reproducibility"))
    st.code("braintrace validate\npython examples/verify_examples.py", language="bash")
    st.write(f"Model lock: `{EXPECTED_MODEL_LOCK_ID}`; locked artifacts: 8/8.")

    st.header(tr("引用、版本与归档", "Citation / version / archive"))
    st.write(
        tr(
            "当前可执行软件版本为 v0.1.16。v0.1.15 是不可变历史归档；v0.1.16 的 release metadata 将在新的不可变归档创建后同步。",
            "The current executable software version is v0.1.16. v0.1.15 remains an immutable historical archive; v0.1.16 release metadata will be synchronized after its immutable archive is created.",
        )
    )
