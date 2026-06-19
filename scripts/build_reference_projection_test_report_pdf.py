#!/usr/bin/env python
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.pdfbase import pdfmetrics
from reportlab.platypus import (
    Image,
    KeepTogether,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)
from reportlab.graphics.shapes import Drawing, String
from reportlab.graphics.charts.barcharts import VerticalBarChart
from reportlab.graphics.charts.legends import Legend


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "bo2023_reference_projection_20260616_cleaned_symbols"
OUTDIR = ROOT / "output" / "pdf"
ASSET_DIR = OUTDIR / "assets"
PDF_PATH = OUTDIR / "reference_projection_test_report_zh.pdf"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def fnum(value: Any, digits: int = 4) -> str:
    if value is None or pd.isna(value):
        return "NA"
    if isinstance(value, (int, float)):
        return f"{float(value):.{digits}f}"
    try:
        return f"{float(value):.{digits}f}"
    except Exception:
        return str(value)


def pct(value: Any) -> str:
    if value is None or pd.isna(value):
        return "NA"
    return f"{float(value) * 100:.1f}%"


def table_data_from_df(df: pd.DataFrame, columns: list[str], headers: list[str]) -> list[list[str]]:
    data = [headers]
    if df.empty:
        return data + [["NA"] * len(headers)]
    for _, row in df.iterrows():
        vals = []
        for col in columns:
            val = row.get(col, "")
            vals.append(fnum(val) if isinstance(val, (int, float)) or str(val).replace(".", "", 1).isdigit() else str(val))
        data.append(vals)
    return data


def make_table(data: list[list[Any]], col_widths: list[float] | None = None, font_size: int = 8) -> Table:
    tbl = Table(data, colWidths=col_widths, repeatRows=1)
    tbl.setStyle(
        TableStyle(
            [
                ("FONTNAME", (0, 0), (-1, -1), "STSong-Light"),
                ("FONTSIZE", (0, 0), (-1, -1), font_size),
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#233142")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#B7C0CA")),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F4F7FA")]),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("ALIGN", (1, 1), (-1, -1), "CENTER"),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    return tbl


def metric_lookup(df: pd.DataFrame, route_col: str, route: str, metric: str) -> float | None:
    if df.empty or route_col not in df.columns or metric not in df.columns:
        return None
    sub = df[df[route_col].astype(str).eq(route)]
    if sub.empty:
        return None
    return float(sub.iloc[0][metric])


def bar_chart(title: str, categories: list[str], series: list[tuple[str, list[float]]], width: int = 470, height: int = 285) -> Drawing:
    drawing = Drawing(width, height)
    chart = VerticalBarChart()
    chart.x = 45
    chart.y = 88
    chart.height = 135
    chart.width = width - 95
    chart.data = [vals for _, vals in series]
    chart.categoryAxis.categoryNames = categories
    chart.categoryAxis.labels.boxAnchor = "ne"
    chart.categoryAxis.labels.angle = 25
    chart.categoryAxis.labels.fontName = "Helvetica"
    chart.categoryAxis.labels.fontSize = 7
    chart.valueAxis.valueMin = 0
    chart.valueAxis.valueMax = 1.0
    chart.valueAxis.valueStep = 0.2
    chart.valueAxis.labels.fontName = "Helvetica"
    chart.valueAxis.labels.fontSize = 7
    palette = ["#2F6F9F", "#D0833D", "#4C956C", "#8D5A97"]
    for idx in range(len(series)):
        chart.bars[idx].fillColor = colors.HexColor(palette[idx % len(palette)])
    drawing.add(chart)
    drawing.add(String(width / 2, height - 20, title, textAnchor="middle", fontName="Helvetica-Bold", fontSize=11))
    legend = Legend()
    legend.x = 55
    legend.y = 38
    legend.fontName = "Helvetica"
    legend.fontSize = 7
    legend.colorNamePairs = [
        (colors.HexColor(palette[idx % len(palette)]), name) for idx, (name, _) in enumerate(series)
    ]
    drawing.add(legend)
    return drawing


def header_footer(canvas, doc) -> None:
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(colors.HexColor("#5A6470"))
    canvas.drawString(1.4 * cm, 0.9 * cm, "Reference-projected VSD validation report")
    canvas.drawRightString(A4[0] - 1.4 * cm, 0.9 * cm, f"Page {doc.page}")
    canvas.restoreState()


def img(path: Path, width: float = 15.5 * cm) -> Image | None:
    if not path.exists():
        return None
    image = Image(str(path))
    ratio = image.imageHeight / float(image.imageWidth)
    image.drawWidth = width
    image.drawHeight = width * ratio
    return image


def build_pdf() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    pdfmetrics.registerFont(UnicodeCIDFont("STSong-Light"))

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="CNTitle", parent=styles["Title"], fontName="STSong-Light", fontSize=22, leading=28, alignment=TA_CENTER, textColor=colors.HexColor("#1E2A36")))
    styles.add(ParagraphStyle(name="CNH1", parent=styles["Heading1"], fontName="STSong-Light", fontSize=16, leading=20, spaceBefore=12, spaceAfter=8, textColor=colors.HexColor("#233142")))
    styles.add(ParagraphStyle(name="CNH2", parent=styles["Heading2"], fontName="STSong-Light", fontSize=12, leading=16, spaceBefore=10, spaceAfter=6, textColor=colors.HexColor("#305674")))
    styles.add(ParagraphStyle(name="CNBody", parent=styles["BodyText"], fontName="STSong-Light", fontSize=9.2, leading=14, alignment=TA_LEFT))
    styles.add(ParagraphStyle(name="Small", parent=styles["BodyText"], fontName="STSong-Light", fontSize=7.5, leading=10, textColor=colors.HexColor("#4E5965")))

    data_audit = read_json(BASE / "data_audit_summary.json")
    projector_qc = read_json(BASE / "projector_qc_summary.json")
    gene_audit = read_json(BASE / "gene_symbol_audit_summary.json")

    loso_network = read_csv(BASE / "bo2023_projected_vsd_loso_route_summary.csv")
    lomo_network = read_csv(BASE / "bo2023_projected_vsd_lomo_route_summary.csv")
    exact_direct_loso = read_csv(BASE / "bo2023_projected_vsd_exact_region_loso_route_summary.csv")
    local_loso = read_csv(BASE / "region_local_rerank_loso_hybrid" / "bo2023_projected_vsd_region_local_rerank_loso_route_summary.csv")
    formal_lomo_net = read_csv(BASE / "formal_three_tier_lomo_hybrid" / "formal_lomo_network_route_metrics.csv")
    formal_lomo_grp = read_csv(BASE / "formal_three_tier_lomo_hybrid" / "formal_lomo_resolution_group_route_metrics.csv")
    formal_lomo_ex = read_csv(BASE / "formal_three_tier_lomo_hybrid" / "formal_lomo_exact_region_route_metrics.csv")
    formal_loso_net = read_csv(BASE / "formal_three_tier_loso_hybrid" / "hybrid_formal_loso_network_route_metrics.csv")
    formal_loso_grp = read_csv(BASE / "formal_three_tier_loso_hybrid" / "hybrid_formal_loso_resolution_group_route_metrics.csv")
    formal_loso_ex = read_csv(BASE / "formal_three_tier_loso_hybrid" / "hybrid_formal_loso_exact_region_route_metrics.csv")
    ahba_direct = read_csv(BASE / "ahba_external_projected_vsd" / "ahba_projected_vsd_external_metrics.csv")
    ahba_formal = read_csv(BASE / "ahba_external_formal_three_tier" / "ahba_formal_three_tier_metrics.csv")
    ahba_special = read_csv(BASE / "ahba_external_formal_three_tier" / "ahba_formal_three_tier_special_labels.csv")
    tcga = read_csv(BASE / "tcga_labeled_hybrid_formal_external" / "tcga_labeled_hybrid_formal_metrics.csv")

    story: list[Any] = []
    story.append(Spacer(1, 2.5 * cm))
    story.append(Paragraph("参数投影器与 Hybrid 三级溯源验证报告", styles["CNTitle"]))
    story.append(Spacer(1, 0.4 * cm))
    story.append(Paragraph("Reference-projected VSD + logCPM local rerank validation summary", styles["CNBody"]))
    story.append(Spacer(1, 0.8 * cm))
    story.append(make_table([
        ["项目", "内容"],
        ["主线路线", "projected VSD Network Top3 beam + logCPM resolution/local exact rerank"],
        ["内部验证", "Bo2023 LOSO / LOMO / exact region / local rerank / formal three-tier"],
        ["外部验证", "AHBA human RNA-seq; TCGA/BraTS MRI-labeled glioma"],
        ["输出目录", str(PDF_PATH.parent)],
    ], [4.0 * cm, 11.5 * cm], 9))
    story.append(Spacer(1, 0.8 * cm))
    story.append(Paragraph("核心结论：projected VSD 单独用于 exact-region 会损失精度；但作为 Network Top3 beam 与 logCPM 的 resolution/local exact rerank 组合后，内部 LOSO exact Top3 达到 0.4533，外部 AHBA exact Top3 达到 0.4286，是当前最稳健的跨域主线候选。", styles["CNBody"]))
    story.append(PageBreak())

    story.append(Paragraph("1. 数据审计与投影器质量", styles["CNH1"]))
    audit_rows = [
        ["指标", "结果"],
        ["Bo2023 common samples", data_audit.get("n_common_samples", "NA")],
        ["Common gene symbols", data_audit.get("n_common_gene_symbols", "NA")],
        ["Locked network genes present", f"{data_audit.get('n_locked_model_genes_present', 'NA')}/{data_audit.get('n_locked_model_genes', 'NA')}"],
        ["Projector median sample Pearson", fnum(projector_qc.get("median_sample_pearson"))],
        ["Projector global Pearson", fnum(projector_qc.get("global_pearson"))],
        ["Projector MAE", fnum(projector_qc.get("mae"))],
        ["Median per-gene R2", fnum(projector_qc.get("median_gene_pearson_r2", projector_qc.get("median_per_gene_r2")))],
        ["Date-like gene symbols after cleaning", gene_audit.get("n_date_like_symbols", "NA")],
    ]
    story.append(make_table(audit_rows, [7.0 * cm, 8.5 * cm], 8.5))
    story.append(Spacer(1, 0.3 * cm))
    qc_img = img(BASE / "figures" / "projector_gene_qc_distributions.png")
    if qc_img:
        story.append(qc_img)
    story.append(Paragraph("解释：Excel 日期样式基因符号已清理；剩余 ENSMFAG fallback 主要是无稳定 human-like symbol 的基因，不影响 locked 200 network gene 覆盖，但外部人类数据 exact 层应避免依赖这些 fallback gene。", styles["Small"]))

    story.append(PageBreak())
    story.append(Paragraph("2. 内部验证总览", styles["CNH1"]))
    story.append(Paragraph("本节按从简单到正式的顺序汇总：Network LOSO/LOMO、direct exact、local rerank、完整三级 hybrid。", styles["CNBody"]))
    story.append(Spacer(1, 0.2 * cm))

    net_rows = [["验证", "route", "Top1", "Top3"]]
    for name, df in [("LOSO network", loso_network), ("LOMO network", lomo_network)]:
        if not df.empty:
            for _, r in df.iterrows():
                net_rows.append([name, r.get("route", ""), fnum(r.get("top1_accuracy", r.get("network_top1"))), fnum(r.get("top3_accuracy", r.get("network_top3")))])
    story.append(make_table(net_rows, [3.4 * cm, 5.2 * cm, 3.1 * cm, 3.1 * cm], 7.8))
    story.append(Spacer(1, 0.3 * cm))
    internal_chart = bar_chart(
        "Internal formal route comparison",
        ["Network Top3", "Group Top3", "Exact Top1", "Exact Top3"],
        [
            (
                "hybrid",
                [
                    metric_lookup(formal_lomo_net, "route_family", "hybrid_projected_network_logcpm_exact", "top3_accuracy") or 0,
                    metric_lookup(formal_lomo_grp, "route_family", "hybrid_projected_network_logcpm_exact", "group_top3_accuracy") or 0,
                    metric_lookup(formal_lomo_ex, "route_family", "hybrid_projected_network_logcpm_exact", "top1_accuracy") or 0,
                    metric_lookup(formal_lomo_ex, "route_family", "hybrid_projected_network_logcpm_exact", "top3_accuracy") or 0,
                ],
            ),
            (
                "logCPM",
                [
                    metric_lookup(formal_lomo_net, "route_family", "logcpm_baseline", "top3_accuracy") or 0,
                    metric_lookup(formal_lomo_grp, "route_family", "logcpm_baseline", "group_top3_accuracy") or 0,
                    metric_lookup(formal_lomo_ex, "route_family", "logcpm_baseline", "top1_accuracy") or 0,
                    metric_lookup(formal_lomo_ex, "route_family", "logcpm_baseline", "top3_accuracy") or 0,
                ],
            ),
            (
                "projected",
                [
                    metric_lookup(formal_lomo_net, "route_family", "projected_vsd", "top3_accuracy") or 0,
                    metric_lookup(formal_lomo_grp, "route_family", "projected_vsd", "group_top3_accuracy") or 0,
                    metric_lookup(formal_lomo_ex, "route_family", "projected_vsd", "top1_accuracy") or 0,
                    metric_lookup(formal_lomo_ex, "route_family", "projected_vsd", "top3_accuracy") or 0,
                ],
            ),
        ],
    )
    story.append(internal_chart)

    story.append(PageBreak())
    story.append(Paragraph("3. 完整三级 Hybrid: LOMO 与 LOSO", styles["CNH1"]))
    formal_rows = [
        ["验证", "Network Top1", "Network Top3", "Group Top1", "Group Top3", "Exact Top1", "Exact Top3"],
        [
            "LOMO hybrid",
            fnum(metric_lookup(formal_lomo_net, "route_family", "hybrid_projected_network_logcpm_exact", "top1_accuracy")),
            fnum(metric_lookup(formal_lomo_net, "route_family", "hybrid_projected_network_logcpm_exact", "top3_accuracy")),
            fnum(metric_lookup(formal_lomo_grp, "route_family", "hybrid_projected_network_logcpm_exact", "group_top1_accuracy")),
            fnum(metric_lookup(formal_lomo_grp, "route_family", "hybrid_projected_network_logcpm_exact", "group_top3_accuracy")),
            fnum(metric_lookup(formal_lomo_ex, "route_family", "hybrid_projected_network_logcpm_exact", "top1_accuracy")),
            fnum(metric_lookup(formal_lomo_ex, "route_family", "hybrid_projected_network_logcpm_exact", "top3_accuracy")),
        ],
        [
            "LOSO hybrid",
            fnum(formal_loso_net.iloc[0].get("top1_accuracy") if not formal_loso_net.empty else None),
            fnum(formal_loso_net.iloc[0].get("top3_accuracy") if not formal_loso_net.empty else None),
            fnum(formal_loso_grp.iloc[0].get("group_top1_accuracy") if not formal_loso_grp.empty else None),
            fnum(formal_loso_grp.iloc[0].get("group_top3_accuracy") if not formal_loso_grp.empty else None),
            fnum(formal_loso_ex.iloc[0].get("top1_accuracy") if not formal_loso_ex.empty else None),
            fnum(formal_loso_ex.iloc[0].get("top3_accuracy") if not formal_loso_ex.empty else None),
        ],
    ]
    story.append(make_table(formal_rows, [2.8 * cm, 2.2 * cm, 2.2 * cm, 2.2 * cm, 2.2 * cm, 2.2 * cm, 2.2 * cm], 7.6))
    story.append(Spacer(1, 0.3 * cm))
    loso_chart = bar_chart(
        "Hybrid formal LOSO vs LOMO",
        ["Network Top3", "Group Top3", "Exact Top1", "Exact Top3"],
        [
            (
                "LOMO",
                [
                    metric_lookup(formal_lomo_net, "route_family", "hybrid_projected_network_logcpm_exact", "top3_accuracy") or 0,
                    metric_lookup(formal_lomo_grp, "route_family", "hybrid_projected_network_logcpm_exact", "group_top3_accuracy") or 0,
                    metric_lookup(formal_lomo_ex, "route_family", "hybrid_projected_network_logcpm_exact", "top1_accuracy") or 0,
                    metric_lookup(formal_lomo_ex, "route_family", "hybrid_projected_network_logcpm_exact", "top3_accuracy") or 0,
                ],
            ),
            (
                "LOSO",
                [
                    formal_loso_net.iloc[0].get("top3_accuracy", 0) if not formal_loso_net.empty else 0,
                    formal_loso_grp.iloc[0].get("group_top3_accuracy", 0) if not formal_loso_grp.empty else 0,
                    formal_loso_ex.iloc[0].get("top1_accuracy", 0) if not formal_loso_ex.empty else 0,
                    formal_loso_ex.iloc[0].get("top3_accuracy", 0) if not formal_loso_ex.empty else 0,
                ],
            ),
        ],
    )
    story.append(loso_chart)
    story.append(Paragraph("完整三级 LOSO 是本线程最新补齐的验证：814 个可评估 fold，Exact Top3=0.4533，Group Top3=0.7236，强于 LOMO 数值；这是当前内部最支持 hybrid 主线的证据。", styles["CNBody"]))

    story.append(PageBreak())
    story.append(Paragraph("4. Exact-region 与 local rerank 证据链", styles["CNH1"]))
    direct_rows = [["验证", "route", "Exact Top1", "Exact Top3", "备注"]]
    for _, r in exact_direct_loso.iterrows() if not exact_direct_loso.empty else []:
        direct_rows.append(["Direct exact LOSO", r.get("route", ""), fnum(r.get("exact_region_top1", r.get("top1_accuracy"))), fnum(r.get("exact_region_top3", r.get("top3_accuracy"))), "locked 200 direct"])
    for _, r in local_loso.iterrows() if not local_loso.empty else []:
        direct_rows.append(["Local rerank LOSO", str(r.get("route", "")).replace("_", " "), fnum(r.get("exact_region_top1")), fnum(r.get("exact_region_top3")), "Top3 network beam + local genes"])
    story.append(make_table(direct_rows, [3.2 * cm, 6.7 * cm, 2.0 * cm, 2.0 * cm, 2.8 * cm], 6.8))
    story.append(Spacer(1, 0.3 * cm))
    story.append(Paragraph("趋势：direct exact 最弱；加入 network beam 和 region-local rerank 后 exact Top3 明显提升；完整三级再加入 resolution group 和 Top50/Top100 fusion 后，hybrid LOSO exact Top3 进一步升到 0.4533。", styles["CNBody"]))

    story.append(PageBreak())
    story.append(Paragraph("5. 外部验证: AHBA 与 TCGA/BraTS", styles["CNH1"]))
    ahba_rows = [["route", "Network Top1", "Network Top3", "Group Top3", "Exact Top1", "Exact Top3"]]
    for _, r in ahba_formal.iterrows() if not ahba_formal.empty else []:
        ahba_rows.append([
            str(r.get("route", "")),
            fnum(r.get("network_top1_accuracy_coarse")),
            fnum(r.get("network_top3_accuracy_coarse")),
            fnum(r.get("group_top3_accuracy_exact_mapped")),
            fnum(r.get("region_top1_accuracy_exact_mapped")),
            fnum(r.get("region_top3_accuracy_exact_mapped")),
        ])
    story.append(Paragraph("AHBA human RNA-seq formal three-tier", styles["CNH2"]))
    story.append(make_table(ahba_rows, [5.8 * cm, 2.0 * cm, 2.0 * cm, 2.0 * cm, 2.0 * cm, 2.0 * cm], 7))
    story.append(Spacer(1, 0.25 * cm))
    ahba_chart = bar_chart(
        "AHBA formal external validation",
        ["Network Top1", "Network Top3", "Group Top3", "Exact Top3"],
        [
            (
                str(r.get("route", ""))[:16],
                [
                    float(r.get("network_top1_accuracy_coarse", 0)),
                    float(r.get("network_top3_accuracy_coarse", 0)),
                    float(r.get("group_top3_accuracy_exact_mapped", 0)),
                    float(r.get("region_top3_accuracy_exact_mapped", 0)),
                ],
            )
            for _, r in ahba_formal.iterrows()
        ] if not ahba_formal.empty else [("NA", [0, 0, 0, 0])],
    )
    story.append(ahba_chart)

    story.append(PageBreak())
    story.append(Paragraph("TCGA/BraTS MRI-labeled glioma", styles["CNH2"]))
    tcga_rows = [["route", "Network Top1", "Network Top3", "Lobe Top1", "Lobe Top3", "Broad Top3"]]
    for _, r in tcga.iterrows() if not tcga.empty else []:
        tcga_rows.append([
            str(r.get("route", "")),
            fnum(r.get("network_top1_accuracy")),
            fnum(r.get("network_top3_accuracy")),
            fnum(r.get("lobe_top1_accuracy")),
            fnum(r.get("lobe_top3_accuracy")),
            fnum(r.get("broad_top3_accuracy")),
        ])
    story.append(make_table(tcga_rows, [5.8 * cm, 2.0 * cm, 2.0 * cm, 2.0 * cm, 2.0 * cm, 2.0 * cm], 7))
    story.append(Spacer(1, 0.25 * cm))
    tcga_img = img(BASE / "tcga_labeled_hybrid_formal_external" / "tcga_labeled_hybrid_formal_accuracy.png")
    if tcga_img:
        story.append(tcga_img)
    story.append(Paragraph("注意：TCGA/BraTS 的 MRI truth 是 human atlas labels，不是 Bo2023 macaque region ID，因此报告 network/lobe/broad accuracy，不报告 Bo2023 exact-region accuracy。", styles["Small"]))

    story.append(PageBreak())
    story.append(Paragraph("6. 结论与下一步", styles["CNH1"]))
    conclusions = [
        ["结论", "证据"],
        ["projected VSD 适合作为 Network Top3 beam", "内部 LOSO Top3=0.9238；AHBA/TCGA 外部保留 network/broad 优势"],
        ["纯 projected VSD 不适合作为 exact 主路线", "formal LOMO pure projected Exact Top3=0.3793，低于 hybrid/logCPM"],
        ["hybrid 是当前跨域主线候选", "内部完整三级 LOSO Exact Top3=0.4533；AHBA Exact Top3=0.4286"],
        ["logCPM 仍是强内部基线", "formal LOMO logCPM Group Top3=0.7057，略高于 hybrid=0.6909"],
        ["外部 exact 需按标签类型解释", "AHBA 可做 mapped exact；TCGA/BraTS 只能做 network/lobe/broad"],
    ]
    story.append(make_table(conclusions, [5.0 * cm, 10.5 * cm], 8.5))
    story.append(Spacer(1, 0.4 * cm))
    story.append(Paragraph("推荐主线：projected VSD Network Top3 beam + logCPM resolution/local exact rerank。短期下一步应固定该路线为外部验证默认路线，并对 AHBA/TCGA 的误差标签做按类别审计；长期应加入人-猕猴 region label 映射层，减少 exact-region 的跨物种标签偏差。", styles["CNBody"]))

    doc = SimpleDocTemplate(
        str(PDF_PATH),
        pagesize=A4,
        rightMargin=1.4 * cm,
        leftMargin=1.4 * cm,
        topMargin=1.5 * cm,
        bottomMargin=1.4 * cm,
        title="Reference-projected VSD validation report",
        author="Codex",
    )
    doc.build(story, onFirstPage=header_footer, onLaterPages=header_footer)


if __name__ == "__main__":
    build_pdf()
    print(PDF_PATH)
