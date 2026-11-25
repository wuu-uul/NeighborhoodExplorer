"""交互式 HTML 基因邻域图生成函数。"""

from pathlib import Path
from typing import Iterable

import plotly.graph_objects as go

from .extractor import GeneFeature


_STRAND_COLORS = {
    "+": "#4C78A8",
    "-": "#F2A541",
    "?": "#8C8C8C",
}
_TARGET_COLOR = "#D62728"


def _arrow_polygon(
    gene: GeneFeature, center_y: float, half_height: float
) -> tuple[list[float], list[float]]:
    """根据链方向构造箭头多边形；未知方向使用矩形表示。"""
    start = float(gene.start)
    end = float(gene.end)
    length = max(end - start, 1.0)

    if gene.strand == "?":
        # 链向未知时不伪造箭头方向，使用普通矩形表示 CDS 区间。
        return (
            [start, end, end, start, start],
            [
                center_y - half_height,
                center_y - half_height,
                center_y + half_height,
                center_y + half_height,
                center_y - half_height,
            ],
        )
    if gene.strand == "+":
        # 正链箭头朝右，箭头头部宽度随基因长度调整并限制在合理范围内。
        head_length = min(max(length * 0.25, 1.0), length * 0.45)
        x_values = [
            start,
            end - head_length,
            end,
            end - head_length,
            start,
            start,
        ]
    elif gene.strand == "-":
        # 负链箭头朝左，与正链使用镜像的多边形顶点顺序。
        head_length = min(max(length * 0.25, 1.0), length * 0.45)
        x_values = [
            start + head_length,
            end,
            end,
            start + head_length,
            start,
            start + head_length,
        ]
    y_values = [
        center_y - half_height,
        center_y - half_height,
        center_y,
        center_y + half_height,
        center_y + half_height,
        center_y - half_height,
    ]
    return x_values, y_values


def write_html_visualization(
    genes: Iterable[GeneFeature],
    output_prefix: str | Path,
    title: str = "基因邻域图",
) -> Path:
    """生成独立可交互的 ``<output_prefix>.html`` 基因邻域图。

    每个 CDS 按其基因组坐标绘制为箭头：正链朝右、负链朝左、链向未知时
    显示为矩形。目标基因使用红色突出；悬停箭头可查看基因名、坐标、链向、
    长度和记录编号。HTML 内嵌 Plotly.js，不需要联网加载绘图库。
    """
    gene_list = list(genes)
    if not gene_list:
        raise ValueError("没有可供绘图的基因数据")

    record_ids = {gene.record_id for gene in gene_list}
    if len(record_ids) > 1:
        raise ValueError("同一张邻域图只能包含一个 GenBank 序列记录")

    output_path = Path(f"{Path(output_prefix)}.html")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    figure = go.Figure()
    center_y = 0.25
    half_height = 0.18
    label_x: list[float] = []
    labels: list[str] = []

    for gene in gene_list:
        x_values, y_values = _arrow_polygon(gene, center_y, half_height)
        color = _TARGET_COLOR if gene.is_target else _STRAND_COLORS.get(
            gene.strand, _STRAND_COLORS["?"]
        )
        gene_name = gene.display_name
        locus_tag = gene.locus_tag or "无"
        gene_qualifier = gene.gene or "无"

        figure.add_trace(
            go.Scatter(
                x=x_values,
                y=y_values,
                mode="lines",
                fill="toself",
                fillcolor=color,
                line={"color": "#333333", "width": 1},
                name=gene_name,
                showlegend=False,
                customdata=[
                    [
                        gene_name,
                        locus_tag,
                        gene_qualifier,
                        gene.start,
                        gene.end,
                        gene.strand,
                        gene.length,
                        gene.record_id,
                    ]
                ]
                * len(x_values),
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "locus_tag: %{customdata[1]}<br>"
                    "gene: %{customdata[2]}<br>"
                    "位置: %{customdata[3]}–%{customdata[4]} bp<br>"
                    "链向: %{customdata[5]}<br>"
                    "长度: %{customdata[6]} bp<br>"
                    "记录: %{customdata[7]}"
                    "<extra></extra>"
                ),
            )
        )
        label_x.append((gene.start + gene.end) / 2)
        labels.append(gene_name)

    # 为图例提供稳定的方向和目标颜色说明；这些空数据 trace 不会遮挡基因箭头。
    for legend_name, legend_color in (
        ("正链 (+)", _STRAND_COLORS["+"]),
        ("负链 (-)", _STRAND_COLORS["-"]),
        ("链向未知 (?)", _STRAND_COLORS["?"]),
        ("目标基因", _TARGET_COLOR),
    ):
        figure.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker={"color": legend_color, "size": 10, "symbol": "square"},
                name=legend_name,
                hoverinfo="skip",
            )
        )

    # 在箭头下方显示基因名，悬停提示则保留完整的结构化属性。
    figure.add_trace(
        go.Scatter(
            x=label_x,
            y=[-0.12] * len(labels),
            mode="text",
            text=labels,
            textposition="bottom center",
            textfont={"size": 10, "color": "#333333"},
            showlegend=False,
            hoverinfo="skip",
        )
    )

    min_start = min(gene.start for gene in gene_list)
    max_end = max(gene.end for gene in gene_list)
    coordinate_span = max(max_end - min_start + 1, 1)
    margin = max(coordinate_span * 0.05, 1)

    figure.update_layout(
        title=title,
        xaxis_title="基因组位置（bp）",
        yaxis={"visible": False, "range": [-0.65, 0.8], "fixedrange": True},
        xaxis={"range": [min_start - margin, max_end + margin], "rangemode": "normal"},
        hovermode="closest",
        template="plotly_white",
        height=360,
        margin={"l": 40, "r": 30, "t": 70, "b": 60},
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.02, "x": 0},
    )

    figure.write_html(
        str(output_path),
        include_plotlyjs=True,
        full_html=True,
        auto_open=False,
    )
    return output_path
