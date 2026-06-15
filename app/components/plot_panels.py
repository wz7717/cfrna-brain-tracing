from __future__ import annotations

import plotly.express as px
import plotly.graph_objects as go

from app.i18n import tr


def make_score_bar(viz):
    fig = px.bar(
        viz,
        x="score",
        y="region",
        orientation="h",
        title=tr("脑区来源综合得分排名", "Integrated source score ranking"),
        color_discrete_sequence=["#2f6df6"],
    )
    fig.update_layout(height=450, margin=dict(l=10, r=10, t=60, b=10))
    return fig


def make_fraction_ci_bar(viz_ci):
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=viz_ci["fraction"],
            y=viz_ci["region"],
            orientation="h",
            error_x=dict(
                type="data",
                symmetric=False,
                array=viz_ci["err_plus"],
                arrayminus=viz_ci["err_minus"],
            ),
            name="fraction",
            marker=dict(color="#74a1ff"),
        )
    )
    fig.update_layout(
        title=tr("去卷积分数与 95% 置信区间", "Deconvolution fraction with 95% confidence interval"),
        height=450,
        margin=dict(l=10, r=10, t=60, b=10),
    )
    return fig


def make_stability_bar(viz_st):
    fig = px.bar(
        viz_st,
        x="stability",
        y="region",
        orientation="h",
        title=tr("Bootstrap 稳定性（Top1 保持频率）", "Bootstrap stability (Top1 retention frequency)"),
        color_discrete_sequence=["#7cd6c1"],
    )
    fig.update_layout(height=450, xaxis=dict(range=[0, 1]), margin=dict(l=10, r=10, t=60, b=10))
    return fig
