"""GenBank CDS 解析与目标基因邻域提取。"""

from dataclasses import dataclass
from pathlib import Path
import warnings

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature


@dataclass(frozen=True)
class GeneFeature:
    """一个 CDS 特征及其报告所需的标准化属性。

    Start 和 End 使用 GenBank 常见的 1-based 闭区间坐标，便于直接阅读；
    Strand 使用 ``+``、``-`` 或 ``?`` 表示正链、负链或未知链方向。
    """

    record_id: str
    locus_tag: str
    gene: str
    start: int
    end: int
    strand: str
    length: int
    is_target: bool = False

    @property
    def display_name(self) -> str:
        """返回适合表格和图形展示的基因名称。"""
        return self.gene or self.locus_tag or "未知基因"


class NeighborhoodExtractor:
    """从 GenBank 文件中定位目标 CDS 并提取其上下游基因。"""

    def __init__(self, input_path: str | Path) -> None:
        self.input_path = Path(input_path)

    @staticmethod
    def _qualifier_value(feature: SeqFeature, key: str) -> str:
        """安全读取 qualifier 的首个非空值；缺失时返回空字符串。"""
        values = feature.qualifiers.get(key, [])
        for value in values:
            normalized = str(value).strip()
            if normalized:
                return normalized
        return ""

    @classmethod
    def _to_gene_feature(
        cls, feature: SeqFeature, record_id: str, is_target: bool = False
    ) -> GeneFeature:
        """将 Biopython 特征转换为稳定、便于输出的数据结构。"""
        location = feature.location
        start = int(location.start) + 1
        end = int(location.end)

        if location.strand == 1:
            strand = "+"
        elif location.strand == -1:
            strand = "-"
        else:
            strand = "?"

        return GeneFeature(
            record_id=record_id,
            locus_tag=cls._qualifier_value(feature, "locus_tag"),
            gene=cls._qualifier_value(feature, "gene"),
            start=start,
            end=end,
            strand=strand,
            length=len(location),
            is_target=is_target,
        )

    def _read_records(self):
        """读取 GenBank 记录，并对文件缺失或空文件给出明确错误。"""
        if not self.input_path.is_file():
            raise FileNotFoundError(f"找不到 GenBank 输入文件：{self.input_path}")

        # 显式管理文件句柄，避免解析器迭代器在文件结束后留下未关闭的句柄。
        with self.input_path.open("r", encoding="utf-8") as input_handle:
            records = list(SeqIO.parse(input_handle, "genbank"))
        if not records:
            raise ValueError(f"GenBank 文件中没有可读取的序列记录：{self.input_path}")
        return records

    def extract(self, target: str, num: int = 5) -> list[GeneFeature]:
        """提取目标基因及其上下游各 ``num`` 个 CDS。

        在多个 GenBank 记录中查找时，每条记录各自排序，邻域不会跨越
        contig 或染色体边界。若目标匹配到多个 CDS，只返回文件中的首个匹配，
        并发出警告。
        """
        normalized_target = target.strip()
        if not normalized_target:
            raise ValueError("目标基因名称不能为空")
        if num < 0:
            raise ValueError("上下游基因数量不能为负数")

        records = self._read_records()
        sorted_cds_by_record: list[tuple[str, list[SeqFeature]]] = []
        matches: list[tuple[int, int]] = []

        # 逐条记录保留 CDS 顺序，避免把不同 contig 的基因邻接在一起。
        for record_index, record in enumerate(records):
            cds_features = [
                feature for feature in record.features if feature.type == "CDS"
            ]
            cds_features.sort(
                key=lambda feature: (
                    int(feature.location.start),
                    int(feature.location.end),
                )
            )
            sorted_cds_by_record.append((record.id, cds_features))

            # 精确匹配 locus_tag 或 gene；按文件记录顺序和基因坐标顺序
            # 搜索，因此第一个匹配项具有确定性。
            for feature_index, feature in enumerate(cds_features):
                if normalized_target in {
                    self._qualifier_value(feature, "locus_tag"),
                    self._qualifier_value(feature, "gene"),
                }:
                    matches.append((record_index, feature_index))

        if not matches:
            raise LookupError(
                f"未在 GenBank 文件中找到目标基因（locus_tag 或 gene）："
                f"{normalized_target}"
            )

        if len(matches) > 1:
            warnings.warn(
                f"目标 {normalized_target!r} 匹配到 {len(matches)} 个 CDS；"
                "将只处理文件中的第一个匹配项。",
                RuntimeWarning,
                stacklevel=2,
            )

        record_index, target_index = matches[0]
        record_id, cds_features = sorted_cds_by_record[record_index]
        first_index = max(0, target_index - num)
        stop_index = min(len(cds_features), target_index + num + 1)

        return [
            self._to_gene_feature(
                feature,
                record_id,
                is_target=(feature_index == target_index),
            )
            for feature_index, feature in enumerate(
                cds_features[first_index:stop_index], start=first_index
            )
        ]
