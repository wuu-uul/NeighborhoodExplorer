"""命令行参数解析模块。"""

import argparse
from pathlib import Path
from typing import Sequence

from .extractor import NeighborhoodExtractor
from .report import write_tsv_report
from .visualization import write_html_visualization


def _non_negative_int(value: str) -> int:
    """将命令行输入转换为非负整数，并在输入无效时给出清晰错误。"""
    try:
        number = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"必须是整数：{value!r}") from exc

    if number < 0:
        raise argparse.ArgumentTypeError("基因邻域数量不能为负数")

    return number


def build_parser() -> argparse.ArgumentParser:
    """创建并配置 NeighborhoodExplorer 的参数解析器。"""
    parser = argparse.ArgumentParser(
        prog="NeighborhoodExplorer",
        description="从 GenBank 文件中提取目标基因及其上下游基因邻域。",
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=Path,
        metavar="GENBANK",
        help="输入的 GenBank 文件路径（必填）",
    )
    parser.add_argument(
        "-t",
        "--target",
        required=True,
        metavar="NAME",
        help="目标基因的 locus_tag 或 gene 名称（必填）",
    )
    parser.add_argument(
        "-n",
        "--num",
        type=_non_negative_int,
        default=5,
        metavar="N",
        help="目标基因上游和下游分别提取的基因数（默认：5）",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("result"),
        metavar="PREFIX",
        help="输出文件统一前缀（默认：result；将生成 .tsv 和 .html 文件）",
    )
    return parser


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """解析命令行参数；argv 参数便于测试或由其他模块调用。"""
    return build_parser().parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """运行完整流程：解析参数、提取邻域并写出 TSV 和 HTML 报告。"""
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        genes = NeighborhoodExtractor(args.input).extract(args.target, args.num)
    except (OSError, ValueError, LookupError) as exc:
        parser.error(str(exc))

    try:
        tsv_path = write_tsv_report(genes, args.output)
        html_path = write_html_visualization(
            genes,
            args.output,
            title=f"目标基因邻域：{args.target}",
        )
    except OSError as exc:
        parser.error(f"写入报告失败：{exc}")

    print(f"邻域基因数：{len(genes)}")
    print(f"TSV 报告：{tsv_path}")
    print(f"HTML 图表：{html_path}")
    return 0
