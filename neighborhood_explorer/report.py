"""结构化报告生成函数。"""

import csv
from pathlib import Path
from typing import Iterable

from .extractor import GeneFeature


def write_tsv_report(
    genes: Iterable[GeneFeature], output_prefix: str | Path
) -> Path:
    """将邻域基因写入 ``<output_prefix>.tsv`` 并返回生成文件路径。

    表格坐标与提取器保持一致，使用 1-based 闭区间。``Gene_Name`` 优先
    写入 gene 名称；若缺失则回退到 locus_tag，再缺失时使用“未知基因”。
    """
    output_path = Path(f"{Path(output_prefix)}.tsv")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # 使用 csv 模块正确处理可能包含制表符、引号等字符的基因名称。
    with output_path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.writer(output_file, delimiter="\t", lineterminator="\n")
        writer.writerow(["Gene_Name", "Start", "End", "Strand", "Length"])
        for gene in genes:
            writer.writerow(
                [
                    gene.display_name,
                    gene.start,
                    gene.end,
                    gene.strand,
                    gene.length,
                ]
            )

    return output_path
