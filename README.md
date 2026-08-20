# PyPLogo

[![PyPI](https://img.shields.io/pypi/v/pyplogo.svg)](https://pypi.org/project/pyplogo/)
[![Python](https://img.shields.io/pypi/pyversions/pyplogo.svg)](https://pypi.org/project/pyplogo/)
[![License: MIT](https://img.shields.io/badge/License-MIT-2f855a.svg)](https://github.com/SLF1303/pyplogo-toolkit/blob/main/LICENSE)

**PyPLogo** 是面向蛋白质二级结构与多序列比对的 Python 可视化工具，提供紧凑、现代、适合论文排版的结构符号，并保留对布局、颜色、图例和生物学标注的细粒度控制。

**PyPLogo** is a Python toolkit for publication-ready visualization of protein secondary structures and multiple-sequence alignments. It combines compact modern structure glyphs with fine control over layout, colors, legends, and biological annotations.

![PyPLogo alignment example using public 1LYZ lysozyme data](https://raw.githubusercontent.com/SLF1303/pyplogo-toolkit/main/docs/images/lysozyme_current_defaults.png)

The public examples below use chain A of hen egg-white lysozyme from [PDB 1LYZ](https://www.rcsb.org/structure/1LYZ). The variant rows and surface-patch markers are illustrative documentation data, not a user dataset or a biological claim.

## 核心能力 / Highlights

- 绘制单条蛋白序列及其 DSSP 二级结构。
- 在多序列比对中选择任意一条序列作为二级结构参考。
- 螺旋使用拼接菱形块，beta 折叠使用扁平箭头，转角使用 `T`、`TT` 或 `TTT`。
- 默认显示主要结构 `H/G/I/E/T`；可选完整模式显示 bend 和 coil。
- 支持二硫键、抗体结合表位、保守性底色、残基编号和多行换行。
- 可调字体、行距、结构间距、颜色、图例、标记、注释轨道和画布尺寸。
- 通过 Matplotlib 导出高分辨率 PNG、矢量 SVG 或 PDF。

- Visualize a single protein sequence with DSSP secondary-structure assignments.
- Select any aligned sequence as the secondary-structure reference.
- Show faceted helix glyphs, flat beta-strand arrows, and compact turn labels.
- Annotate disulfide bonds, antibody epitopes, conservation, and residue numbers.
- Export publication-ready raster or vector figures through Matplotlib.

## 安装 / Installation

```bash
pip install -U pyplogo
```

PyPI 包要求 Python 3.10 或更高版本。使用 PDB/mmCIF 的 DSSP 自动提取功能时，系统还需要可用的 `mkdssp`/DSSP 可执行程序。

The PyPI package requires Python 3.10 or later. Structure extraction from PDB/mmCIF files additionally requires an available `mkdssp`/DSSP executable.

## 单序列快速开始 / Single Sequence

```python
from pyplogo import SecondaryStructureVisualizer

sequence = "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
secondary_structure = "CECCHHHHHHHHHHTTCTTETTECTHHHHHHHHHHHTTESSCEEECTTSCEEETTTTEETTTSCECSSCTTCCCTTCSEGGGGGSSCCHHHHHHHHHHTTSSSGGGGSHHHHHHTTTSCGGGGSTTCCC"

visualizer = SecondaryStructureVisualizer(row_length=50)
figure = visualizer.create_figure(
    sequence=sequence,
    secondary_structure=secondary_structure,
    title="Hen egg-white lysozyme (PDB 1LYZ)",
    show_residue_numbers=True,
)
figure.savefig("secondary_structure.svg", bbox_inches="tight")
```

`sequence` 与 `secondary_structure` 必须等长。PyPLogo 使用 DSSP 字母：`H`、`G`、`I`、`E`、`T`、`S`，coil 使用 `C`。

`sequence` and `secondary_structure` must have the same length. PyPLogo accepts the DSSP codes `H`, `G`, `I`, `E`, `T`, and `S`, with `C` for coil.

## 多序列比对 / Multiple-Sequence Alignment

PyPLogo 接收已经对齐且长度相同的序列；它负责绘图，不负责执行序列比对。gap 可以写作 `-` 或 `.`。

PyPLogo draws an existing equal-length alignment; it does not run the alignment algorithm. Gaps may be written as `-` or `.`.

```python
from pyplogo import SecondaryStructureVisualizer

alignment = {
    "Lysozyme_1LYZ": "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL",
    "Lysozyme_variant_A": "KVFGRCELAVAMKRHGLDNYRGYSLGNWVCAAKFESNSNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL",
    "Lysozyme_variant_B": "KVFGRCELAAAMKRHGQDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL",
}

visualizer = SecondaryStructureVisualizer(
    colors={
        "H": "#39A96B",
        "G": "#D982AA",
        "I": "#6878C9",
        "E": "#D7A354",
        "T": "#8F7CC3",
        "structure_label": "#9A642D",
        "background": "#FFFFFF",
    }
)

figure = visualizer.create_alignment_figure(
    alignment=alignment,
    structure_sequence="Lysozyme_1LYZ",
    # DSSP assignment for the public 1LYZ chain A sequence.
    secondary_structure="CECCHHHHHHHHHHTTCTTETTECTHHHHHHHHHHHTTESSCEEECTTSCEEETTTTEETTTSCECSSCTTCCCTTCSEGGGGGSSCCHHHHHHHHHHTTSSSGGGGSHHHHHHTTTSCGGGGSTTCCC",
    structure_label="1LYZ structure",
    amino_acids_per_line=65,
    structure_mode="major",
    sequence_colors={
        "Lysozyme_1LYZ": "#986B36",
        "Lysozyme_variant_A": "#5B6B7A",
        "Lysozyme_variant_B": "#6F68A8",
    },
    # Disulfides reported for the public 1LYZ chain A structure.
    disulfide_bonds=[(6, 127), (30, 115), (64, 80), (76, 94)],
    epitope_annotations={
        "Illustrative surface patch A": [22, 23, 24, 25],
        "Illustrative surface patch B": [52, 53, 54, 55],
        "Illustrative surface patch C": [101, 102, 103, 104],
    },
    epitope_colors={
        "Illustrative surface patch A": "#8CB9DC",
        "Illustrative surface patch B": "#D59A9A",
        "Illustrative surface patch C": "#A8C565",
    },
    show_legend=True,
    legend_position="bottom",
)
figure.savefig("alignment.svg", bbox_inches="tight")
figure.savefig("alignment.png", dpi=600, bbox_inches="tight")
```

`secondary_structure` 可以使用所选参考序列的**无 gap 长度**，也可以使用完整比对长度。二硫键和表位位置均使用所选参考序列的 **1-based 无 gap 残基编号**。

示例中的三个 surface patch 只是演示注释轨道的占位数据；实际项目中可以替换为实验或结构分析得到的残基位置。

`secondary_structure` may follow either the selected sequence's **ungapped length** or the full alignment length. Disulfide and epitope positions use **1-based ungapped residue coordinates** from the selected reference sequence.

## 结构符号 / Structure Glyphs

| Code | DSSP structure | PyPLogo representation |
| --- | --- | --- |
| `H` | alpha-helix | 绿色拼接菱形块 / green interlocking facets |
| `G` | 3_10-helix | 粉色拼接菱形块 / pink interlocking facets |
| `I` | pi-helix | 蓝紫色拼接菱形块 / blue-violet interlocking facets |
| `E` | beta-strand | 扁平箭头 / flat arrow |
| `T` | hydrogen-bonded turn | `T`、`TT`、`TTT` |
| `S` | bend | 折线，仅完整模式显示 / bend line, full mode only |
| `C` | coil | 细横线，仅完整模式显示 / thin line, full mode only |

`T`、`TT` 和 `TTT` 对应连续 1、2、3 个或更多 turn 残基。默认 `structure_mode="major"` 显示 `H/G/I/E/T`；设置 `structure_mode="all"` 后额外显示 `S/C` 和 gap 连接符。

`T`, `TT`, and `TTT` represent runs of one, two, or at least three turn residues. The default `structure_mode="major"` shows `H/G/I/E/T`; `structure_mode="all"` additionally shows bends, coils, and gap connectors.

## 自定义参数 / Customization

| Category | Frequently used parameters |
| --- | --- |
| Reference and wrapping | `structure_sequence`, `structure_label`, `amino_acids_per_line`, `figsize` |
| Typography and spacing | `sequence_font_size`, `label_font_size`, `sequence_pitch`, `structure_offset` |
| Structure display | `structure_mode`, `structure_gap`, constructor `colors` (`H/G/I/E/T/S/C`, labels, background) |
| Conservation | `show_conservation`, `conservation_threshold`, `conservation_full_color`, `conservation_partial_color` |
| Residue scaffolding | `show_residue_numbers`, `sequence_colors` |
| Structure legend | `show_legend`, `legend_position`, `legend_columns`, `legend_font_size` |
| Disulfides | `disulfide_bonds`, `show_disulfide_bonds`, colors, line length/width, label and line offsets |
| Epitope tracks | `epitope_annotations`, `epitope_colors`, marker type/size, legend position/columns |
| Annotation layout | `annotation_offset`, `annotation_track_spacing` |

The defaults provide the base publication style; every group above can be overridden per figure without changing the renderer.

![PyPLogo customization gallery using public lysozyme data](https://raw.githubusercontent.com/SLF1303/pyplogo-toolkit/main/docs/images/lysozyme_customization_gallery.png)

## PDB/mmCIF 提取 / Structure Extraction

```python
from pyplogo.extractors.structure_based import StructureExtractor
from pyplogo import SecondaryStructureVisualizer

data = StructureExtractor().from_pdb("1LYZ.pdb", chain_id="A")

visualizer = SecondaryStructureVisualizer()
figure = visualizer.create_figure(
    sequence=data.sequence,
    secondary_structure=data.secondary_structure,
    disulfide_bonds=data.disulfide_bonds,
    title="Hen egg-white lysozyme (PDB 1LYZ), chain A",
)
figure.savefig("chain_A.pdf", bbox_inches="tight")
```

`from_cif("protein.cif", chain_id="A")` 提供相同的 mmCIF 工作流。

## Links

- GitHub: https://github.com/SLF1303/pyplogo-toolkit
- PyPI: https://pypi.org/project/pyplogo/
- Issues: https://github.com/SLF1303/pyplogo-toolkit/issues

PyPLogo is released under the [MIT License](https://github.com/SLF1303/pyplogo-toolkit/blob/main/LICENSE).
