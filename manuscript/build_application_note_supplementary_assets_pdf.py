from __future__ import annotations

from pathlib import Path

import pandas as pd
from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    Image,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "manuscript" / "figures_application_note_20260614" / "supplementary"
TABLE_DIR = ROOT / "manuscript" / "tables_application_note_20260614" / "supplementary"
OUTPUT = ROOT / "manuscript" / "Bioinformatics_Application_Note_Supplementary_Figures_Tables_20260614.pdf"

PAGE = landscape(A4)
PAGE_W, PAGE_H = PAGE
MARGIN = 12 * mm
CONTENT_W = PAGE_W - 2 * MARGIN
CONTENT_H = PAGE_H - 2 * MARGIN

pdfmetrics.registerFont(TTFont("SimHei", r"C:\Windows\Fonts\simhei.ttf"))

styles = getSampleStyleSheet()
title_style = ParagraphStyle(
    "TitleCN",
    parent=styles["Title"],
    fontName="SimHei",
    fontSize=19,
    leading=24,
    alignment=TA_CENTER,
    textColor=colors.HexColor("#17365D"),
    spaceAfter=10,
)
subtitle_style = ParagraphStyle(
    "SubtitleCN",
    parent=styles["BodyText"],
    fontName="SimHei",
    fontSize=11,
    leading=16,
    alignment=TA_CENTER,
    textColor=colors.HexColor("#4B5563"),
)
heading_style = ParagraphStyle(
    "HeadingCN",
    parent=styles["Heading1"],
    fontName="SimHei",
    fontSize=13,
    leading=16,
    textColor=colors.HexColor("#1F4E79"),
    spaceAfter=5,
)
caption_style = ParagraphStyle(
    "Caption",
    parent=styles["BodyText"],
    fontName="Helvetica",
    fontSize=8.2,
    leading=10.5,
    textColor=colors.HexColor("#374151"),
    spaceBefore=4,
)
cell_style = ParagraphStyle(
    "Cell",
    parent=styles["BodyText"],
    fontName="Helvetica",
    fontSize=7.2,
    leading=8.8,
)
header_style = ParagraphStyle(
    "HeaderCell",
    parent=cell_style,
    fontName="Helvetica-Bold",
    textColor=colors.white,
    alignment=TA_CENTER,
)


FIGURES = [
    (
        "Figure S1. Internal validation across anatomical resolutions",
        "FigureS1_internal_validation.png",
        "LOSO and leave-one-monkey-out performance at Network, Region Group and Exact Region levels.",
    ),
    (
        "Figure S2. Pairwise-rescue threshold and donor-isolated validation",
        "FigureS2_threshold_lomo.png",
        "Retrospective margin-threshold audit and leave-one-monkey-out generalization.",
    ),
    (
        "Figure S3. AHBA cross-species validation",
        "FigureS3_AHBA_external_validation.png",
        "Harmonized coarse-label performance and the limitations of exact-region transfer.",
    ),
    (
        "Figure S4. Glioma and liquid-biopsy transfer diagnostics",
        "FigureS4_glioma_liquid_biopsy.png",
        "Paired glioma transcriptome-MRI results and liquid-biopsy transfer limits. Biofluid cohorts lacked imaging truth.",
    ),
    (
        "Figure S5. Exploratory adaptation sensitivity",
        "FigureS5_adaptation_sensitivity.png",
        "Endpoint-dependent adaptation results. Adapted routes are supplementary and are not part of the locked workflow.",
    ),
]

TABLES = [
    (
        "Table S1. Datasets and analysis roles",
        "TableS1_datasets_and_roles.csv",
        [31, 45, 32, 54, 57],
    ),
    (
        "Table S2. Internal validation metrics",
        "TableS2_internal_validation.csv",
        [20, 27, 25, 25, 39, 25, 39],
    ),
    (
        "Table S3. External validation metrics",
        "TableS3_external_validation.csv",
        [30, 58, 20, 42, 76],
    ),
    (
        "Table S4. Domain-adaptation sensitivity",
        "TableS4_domain_adaptation_sensitivity.csv",
        [37, 73, 42, 38, 31, 31, 31],
    ),
    (
        "Table S5. Liquid-biopsy transfer stress tests",
        "TableS5_liquid_biopsy_stress_tests.csv",
        [28, 43, 39, 49, 47, 79],
    ),
]


def page_number(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(colors.HexColor("#6B7280"))
    canvas.drawRightString(PAGE_W - MARGIN, 6 * mm, f"Page {doc.page}")
    canvas.restoreState()


def fitted_image(path: Path, max_w: float, max_h: float) -> Image:
    with PILImage.open(path) as img:
        width, height = img.size
    scale = min(max_w / width, max_h / height)
    return Image(str(path), width=width * scale, height=height * scale)


def paragraph_cell(value, style=cell_style):
    text = str(value).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return Paragraph(text, style)


def build_table(path: Path, widths_mm: list[float]) -> Table:
    frame = pd.read_csv(path, dtype=str).fillna("")
    data = [[paragraph_cell(column, header_style) for column in frame.columns]]
    for row in frame.itertuples(index=False, name=None):
        data.append([paragraph_cell(value) for value in row])
    total = sum(widths_mm)
    widths = [CONTENT_W * value / total for value in widths_mm]
    table = Table(data, colWidths=widths, repeatRows=1, hAlign="CENTER")
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1F4E79")),
                ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#9CA3AF")),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F3F6F9")]),
            ]
        )
    )
    return table


def main() -> None:
    doc = SimpleDocTemplate(
        str(OUTPUT),
        pagesize=PAGE,
        leftMargin=MARGIN,
        rightMargin=MARGIN,
        topMargin=MARGIN,
        bottomMargin=MARGIN,
        title="Bioinformatics Application Note Supplementary Figures and Tables",
        author="",
    )
    story = [
        Spacer(1, 35 * mm),
        Paragraph("Bioinformatics Application Note补充图表", title_style),
        Paragraph("cfRNA-BrainTrace Supplementary Figures and Tables", subtitle_style),
        Spacer(1, 10 * mm),
        Paragraph("包含Figure S1-S5及Table S1-S5", subtitle_style),
        Paragraph("生成日期：2026年6月14日", subtitle_style),
        PageBreak(),
    ]

    for index, (title, filename, caption) in enumerate(FIGURES):
        story.append(Paragraph(title, heading_style))
        story.append(Spacer(1, 2 * mm))
        story.append(fitted_image(FIG_DIR / filename, CONTENT_W, CONTENT_H - 29 * mm))
        story.append(Paragraph(caption, caption_style))
        story.append(PageBreak())

    for index, (title, filename, widths) in enumerate(TABLES):
        story.append(Paragraph(title, heading_style))
        story.append(Spacer(1, 3 * mm))
        story.append(build_table(TABLE_DIR / filename, widths))
        if index < len(TABLES) - 1:
            story.append(PageBreak())

    doc.build(story, onFirstPage=page_number, onLaterPages=page_number)
    print(OUTPUT)


if __name__ == "__main__":
    main()
