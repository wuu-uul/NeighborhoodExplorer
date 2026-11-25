# NeighborhoodExplorer

[简体中文](README.md) | [English](README.en.md)

NeighborhoodExplorer 是一个 Python 命令行工具，可从 GenBank 注释文件中按
`locus_tag` 或 `gene` 名称查找目标 CDS，提取其上下游基因邻域，并生成 TSV
数据表和可交互的 HTML 基因结构图。

## 许可证

本项目使用 [MIT 许可证](LICENSE)。

## 功能

- 使用 Biopython 解析 GenBank 文件中的 CDS 特征。
- 按 CDS 起始位置排序，并在目标所在的序列记录内提取邻域。
- 目标基因上游、下游各提取指定数量的 CDS；遇到记录边界时只返回可用基因。
- 支持通过 `locus_tag` 或 `gene` 精确匹配目标。多个 CDS 匹配时发出警告并使用第一个匹配项。
- 输出包含基因名称、坐标、链方向和长度的 TSV 文件。
- 输出可离线打开的 Plotly HTML 图：正链箭头朝右、负链箭头朝左，目标 CDS 高亮，悬停可查看详细信息。

## 环境要求

- Python 3.10 或更高版本
- Biopython
- Plotly

## 安装

在项目根目录执行：

```powershell
python -m pip install -r requirements.txt
```

## 快速开始

```powershell
python -m gene --input genome.gbk --target ABC_0001 --num 5 --output results/neighborhood
```

也可以使用短参数：

```powershell
python -m gene -i genome.gbk -t ABC_0001 -n 5 -o results/neighborhood
```

运行成功后会生成：

```text
results/
├── neighborhood.tsv
└── neighborhood.html
```

`--output` 是统一文件前缀，不需要包含 `.tsv` 或 `.html` 扩展名。未指定时，
前缀默认为当前目录下的 `result`。

## 命令行参数

| 参数 | 长参数 | 必填 | 默认值 | 说明 |
|---|---|---:|---|---|
| `-i` | `--input` | 是 | — | 输入 GenBank 文件路径 |
| `-t` | `--target` | 是 | — | 目标基因的 `locus_tag` 或 `gene` 名称 |
| `-n` | `--num` | 否 | `5` | 目标 CDS 上游和下游分别提取的 CDS 数量；必须是非负整数 |
| `-o` | `--output` | 否 | `result` | 输出文件统一前缀 |

例如，只提取目标 CDS 本身：

```powershell
python -m gene -i genome.gbk -t ABC_0001 -n 0 -o results/target_only
```

查看完整参数说明：

```powershell
python -m gene --help
```

## 输出说明

### TSV 数据表

TSV 文件使用制表符分隔，列顺序如下：

| 列名 | 说明 |
|---|---|
| `Gene_Name` | `gene` 名称；缺失时使用 `locus_tag`，两者都缺失时显示“未知基因” |
| `Start` | 基因起始位置，1-based 坐标 |
| `End` | 基因结束位置，1-based 闭区间 |
| `Strand` | `+` 表示正链，`-` 表示负链，`?` 表示方向未知 |
| `Length` | CDS 特征长度，单位为 bp |

### HTML 邻域图

HTML 文件包含 Plotly.js，可以在浏览器中直接打开，无需联网加载绘图库。
每个 CDS 按其坐标绘制；悬停可查看名称、`locus_tag`、坐标、链方向、长度和
GenBank 记录编号。目标基因以红色高亮，图例说明链方向颜色。

## 邻域与匹配规则

- 只将 `type="CDS"` 的注释特征视为基因。
- 邻域按基因组坐标排序，不按链方向反转顺序。
- 多条记录（例如多个 contig）分别处理，邻域不会跨记录边界。
- 目标名称与 `locus_tag` 或 `gene` 的值进行精确、区分大小写的匹配。
- 如果有多个 CDS 匹配，程序会发出警告，并按 GenBank 记录顺序及坐标顺序选择第一个。
- `--num N` 表示目标两侧各最多取 N 个 CDS，因此完整邻域最多包含 `2N + 1` 个 CDS。
- TSV 和图表中的起止位置使用便于阅读的 1-based 坐标；GenBank 特征的结束位置按闭区间展示。

## 项目结构

```text
gene/
├── __init__.py       # 包信息
├── __main__.py       # python -m gene 入口
├── cli.py            # 命令行参数与完整执行流程
├── extractor.py      # GenBank 解析及邻域提取
├── report.py         # TSV 报告生成
└── visualization.py  # HTML 基因邻域图生成
requirements.txt      # Python 运行依赖
```

## 常见问题

### 找不到目标基因

确认目标值与 GenBank CDS 的 `locus_tag` 或 `gene` 注释完全一致。匹配区分大小写，
并且只搜索 CDS 特征。

### 输出基因数量少于预期

目标靠近序列记录开头或结尾时，程序只会提取实际存在的邻域基因；也不会跨不同
GenBank 记录提取基因。

### GenBank 文件或报告写入错误

检查输入路径和文件权限。输入文件必须包含可解析的 GenBank 序列记录；输出目录
会自动创建，但输出位置仍需具有写入权限。

## Python API 示例

除了命令行，也可以在 Python 代码中调用提取器和报告函数：

```python
from gene.extractor import NeighborhoodExtractor
from gene.report import write_tsv_report
from gene.visualization import write_html_visualization

genes = NeighborhoodExtractor("genome.gbk").extract("ABC_0001", num=5)
write_tsv_report(genes, "results/neighborhood")
write_html_visualization(genes, "results/neighborhood")
```
