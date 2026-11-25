# NeighborhoodExplorer

[简体中文](README.md) | [English](README.en.md)

NeighborhoodExplorer is a Python command-line tool that finds a target CDS by
`locus_tag` or `gene` in an annotated GenBank file, extracts its neighboring
genes, and generates a TSV table and an interactive HTML gene map.

## License

This project is licensed under the [MIT License](LICENSE).

## Features

- Parses CDS features in GenBank files with Biopython.
- Sorts CDS features by genomic start coordinate and extracts neighbors within
  the target's sequence record.
- Extracts up to the requested number of CDS features on each side of the
  target; fewer are returned when a record boundary is reached.
- Matches the target exactly against `locus_tag` or `gene`. If multiple CDS
  features match, a warning is issued and the first match is used.
- Writes a TSV report with gene names, coordinates, strand, and length.
- Creates a standalone Plotly HTML map: positive-strand genes point right,
  negative-strand genes point left, and the target CDS is highlighted.
  Hover over a gene to see its details.

## Requirements

- Python 3.10 or later
- Biopython
- Plotly

## Installation

Run this command from the project directory:

```powershell
python -m pip install -r requirements.txt
```

## Quick start

New to PowerShell? See the [Windows step-by-step guide in Chinese](docs/使用指南.md).

```powershell
python -m neighborhood_explorer --input genome.gbk --target ABC_0001 --num 5 --output results/neighborhood
```

Short options are also available:

```powershell
python -m neighborhood_explorer -i genome.gbk -t ABC_0001 -n 5 -o results/neighborhood
```

On success, the tool creates:

```text
results/
├── neighborhood.tsv
└── neighborhood.html
```

`--output` is a shared file prefix; do not include `.tsv` or `.html`. If omitted,
the prefix defaults to `result` in the current directory.

## Command-line options

| Short option | Long option | Required | Default | Description |
|---|---|---:|---|---|
| `-i` | `--input` | Yes | — | Path to the input GenBank file |
| `-t` | `--target` | Yes | — | Target gene's `locus_tag` or `gene` name |
| `-n` | `--num` | No | `5` | Number of CDS features to extract upstream and downstream; must be a non-negative integer |
| `-o` | `--output` | No | `result` | Shared output file prefix |

To extract only the target CDS, set the flank size to zero:

```powershell
python -m neighborhood_explorer -i genome.gbk -t ABC_0001 -n 0 -o results/target_only
```

For all available options:

```powershell
python -m neighborhood_explorer --help
```

## Output files

### TSV table

The TSV file is tab-delimited and uses the following column order:

| Column | Description |
|---|---|
| `Gene_Name` | The `gene` name; falls back to `locus_tag`, then `未知基因` if both are missing |
| `Start` | Gene start coordinate, using 1-based coordinates |
| `End` | Gene end coordinate, using a 1-based closed interval |
| `Strand` | `+` for positive, `-` for negative, and `?` for unknown strand |
| `Length` | CDS feature length in bp |

### HTML neighborhood map

The HTML file bundles Plotly.js and can be opened in a browser without
downloading the plotting library at runtime. CDS features are drawn at their
genomic coordinates. Hover over an arrow to see the gene name, `locus_tag`,
coordinates, strand, length, and GenBank record ID. The target gene is
highlighted in red, and the legend explains strand colors.

## Matching and neighborhood rules

- Only annotated features with `type="CDS"` are treated as genes.
- Features are ordered by genomic coordinates; order is not reversed for genes
  on the negative strand.
- GenBank records (for example, separate contigs) are handled independently.
  Neighborhood extraction does not cross record boundaries.
- The target is matched exactly and case-sensitively against a CDS feature's
  `locus_tag` or `gene`.
- If multiple CDS features match, the tool warns and selects the first one in
  GenBank record order and coordinate order.
- `--num N` means up to N CDS features on each side, so a complete neighborhood
  contains at most `2N + 1` CDS features.
- Start and end coordinates in the TSV and chart are displayed as 1-based
  coordinates, with the end coordinate included.

## Project structure

```text
neighborhood_explorer/
├── __init__.py       # Package metadata
├── __main__.py       # Entry point for python -m neighborhood_explorer
├── cli.py            # Command-line arguments and application workflow
├── extractor.py      # GenBank parsing and neighborhood extraction
├── report.py         # TSV report generation
└── visualization.py  # HTML gene map generation
requirements.txt      # Python runtime dependencies
```

## Troubleshooting

### Target gene not found

Check that the target exactly matches the `locus_tag` or `gene` qualifier of a
CDS feature. Matching is case-sensitive, and only CDS features are searched.

### Fewer genes than expected

If the target is close to the beginning or end of a sequence record, only the
available neighboring CDS features are returned. The tool does not continue
into another GenBank record.

### GenBank input or report-writing errors

Check the input path and file permissions. The input must contain at least one
parseable GenBank sequence record. Output directories are created
automatically, but the process must have permission to write to the output
location.

## Python API

The extractor and report functions can also be used from Python:

```python
from neighborhood_explorer.extractor import NeighborhoodExtractor
from neighborhood_explorer.report import write_tsv_report
from neighborhood_explorer.visualization import write_html_visualization

genes = NeighborhoodExtractor("genome.gbk").extract("ABC_0001", num=5)
write_tsv_report(genes, "results/neighborhood")
write_html_visualization(genes, "results/neighborhood")
```
