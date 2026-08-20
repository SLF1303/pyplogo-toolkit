import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
from matplotlib import gridspec
from typing import Optional, List, Tuple, Dict, Union, Callable, Mapping, Sequence
import warnings
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Polypeptide import is_aa

class SecondaryStructureVisualizer:
    """Create publication-quality secondary structure visualizations with professional scientific style"""
    
    # 默认颜色方案
    DEFAULT_COLORS = {
        'H': '#4A9D5B',      # α-helix - dark green
        'G': '#C97BA0',      # 3₁₀-helix - pink
        'I': '#4169E1',      # π-helix - royal blue
        'E': '#F5D598',      # β-strand - light yellow
        'T': '#C8B6E2',      # turn - light purple
        'S': '#8CD0E3',      # bend - light blue
        'C': '#A69F98',      # coil - gray
        'sequence': '#2D3748', # amino acid sequence
        'structure_label': '#616971', # secondary-structure track label
        'residue_number': '#718096',   # residue index
        'edge': '#4A4643',    # structure edges
        'background': '#FFFAF5', # background color
        'highlight_aa': '#FFE0B2',  # 氨基酸高亮背景
        'highlight_ss': '#E1F5FE',  # 二级结构高亮背景
        'highlight_text': '#D32F2F', # 高亮区域标题颜色
        'legend_text': '#5D534A',   # 图例文本颜色
        'disulfide_text': "#4DFF00"   # 二硫键文本颜色
    }
    
    def __init__(self, 
                 row_length=50,
                 amino_acid_spacing=28,
                 helix_intra_unit_overlap=12,
                 helix_inter_unit_gap=0,
                 parallelogram_width=15,
                 helix_unit_height=8.5,
                 helix_skew=12,
                 beta_body_width_ratio=0.3,
                 beta_head_width_ratio=0.55,
                 beta_head_length_ratio=1.0,
                 aa_font_size=15,
                 aa_font_family='Times New Roman',
                 # 新增比例控制参数
                 fig_width_scale=1.0,
                 fig_height_scale=1.0,
                 aspect_ratio=None,
                 top_margin=0.1,
                 bottom_margin=0.1,
                 left_margin=0.1,
                 right_margin=0.1,
                 horizontal_spacing=0.05,
                 vertical_spacing=0.05,
                 # 新增颜色方案参数
                 colors: Optional[Dict[str, str]] = None,
                 disulfide_font_size: int = 12,
                 disulfide_height_offset: float = 17.0,
                 residue_number_font_size: Optional[int] = None,
                 # 新增图例位置参数
                 legend_position: str = 'right',  # 'bottom' 或 'right'
                 residue_number_position: float = 0.5):  # 0.0=最低, 1.0=最高
                
        """
        Initialize visualizer with professional scientific style
        
        Args:
            row_length: Number of amino acids per row
            amino_acid_spacing: Spacing between amino acids
            helix_intra_unit_overlap: Overlap within helix unit
            helix_inter_unit_gap: Gap between helix units
            parallelogram_width: Width of parallelogram in helix
            helix_unit_height: Height of helix unit
            helix_skew: Skew degree of helix
            beta_body_width_ratio: Beta strand body width ratio
            beta_head_width_ratio: Beta strand head width ratio
            beta_head_length_ratio: Beta strand head length ratio
            aa_font_size: Amino acid font size
            aa_font_family: Amino acid font family
            fig_width_scale: 图形宽度缩放因子 (默认: 1.0)
            fig_height_scale: 图形高度缩放因子 (默认: 1.0)
            aspect_ratio: 宽高比 (None=自动, 或指定如 16/9)
            top_margin: 顶部边距比例 (0-1)
            bottom_margin: 底部边距比例 (0-1)
            left_margin: 左侧边距比例 (0-1)
            right_margin: 右侧边距比例 (0-1)
            horizontal_spacing: 水平间距比例 (0-1)
            vertical_spacing: 垂直间距比例 (0-1)
            colors: 自定义颜色方案字典，覆盖默认颜色
            disulfide_font_size: 二硫键标记字体大小
             disulfide_height_offset: 二硫键标记高度偏移量 (默认: 17.0，值越大位置越低)
             residue_number_font_size: 残基编号字体大小 (默认: None，自动使用 aa_font_size - 2)
             legend_position: 图例位置 ('bottom' 或 'right')
        """
        # 存储所有参数
        self.row_length = row_length
        self.aa_spacing = amino_acid_spacing
        self.aa_font_size = aa_font_size
        self.aa_font_family = aa_font_family
        self.legend_position = legend_position
        self.disulfide_font_size = disulfide_font_size
        self.disulfide_height_offset = disulfide_height_offset
        self.residue_number_font_size = residue_number_font_size if residue_number_font_size is not None else aa_font_size - 2
        self.residue_number_position = residue_number_position
        
        # 螺旋参数
        self.helix_intra_unit_overlap = helix_intra_unit_overlap
        self.helix_inter_unit_gap = helix_inter_unit_gap
        self.parallelogram_width = parallelogram_width
        self.helix_unit_height = helix_unit_height
        self.helix_skew = helix_skew
        
        # Beta折叠参数
        self.beta_body_width_ratio = beta_body_width_ratio
        self.beta_head_width_ratio = beta_head_width_ratio
        self.beta_head_length_ratio = beta_head_length_ratio
        
        # 新增比例参数
        self.fig_width_scale = fig_width_scale
        self.fig_height_scale = fig_height_scale
        self.aspect_ratio = aspect_ratio
        self.top_margin = top_margin
        self.bottom_margin = bottom_margin
        self.left_margin = left_margin
        self.right_margin = right_margin
        self.horizontal_spacing = horizontal_spacing
        self.vertical_spacing = vertical_spacing
        
        # 计算派生参数
        self.helix_unit_total_width = (2 * self.parallelogram_width - 
                                      self.helix_intra_unit_overlap)
        self.position_compensation = (self.helix_unit_total_width - 
                                      self.aa_spacing + 
                                      self.helix_inter_unit_gap) / 2
        
        # 布局参数
        self.row_height = 60 + (aa_font_size - 10) * 2 + max(0, self.helix_unit_height - 8) * 1.5
        self.base_fig_width = 10 + (row_length * self.aa_spacing) / 20
        self.base_fig_height_per_row = 3.0 + (aa_font_size - 10) * 0.15 + max(0, self.helix_unit_height - 8) * 0.08
        
        # 应用比例缩放
        self.base_fig_width *= self.fig_width_scale
        self.base_fig_height_per_row *= self.fig_height_scale
        
        # 结构元素扩展
        self.max_element_extension = {
            'left': 30 + self.helix_skew,
            'right': 40 + self.helix_skew + int(self.aa_spacing * (self.beta_head_length_ratio - 0.5)),
            'top': 40 + (aa_font_size - 10) + self.helix_unit_height,
            'bottom': 35 + (aa_font_size - 10)
        }
        
        # Beta折叠尺寸
        self.beta_head_width = self.aa_spacing * self.beta_head_width_ratio
        self.beta_head_std_length = self.aa_spacing * self.beta_head_length_ratio
        self.beta_body_width = self.aa_spacing * self.beta_body_width_ratio
        self.beta_body_min_length = self.aa_spacing * 0.7
        
        self.smooth_points = 150
        
        # 颜色方案 - 使用自定义颜色或默认颜色
        self.colors = self.DEFAULT_COLORS.copy()
        if colors:
            # 只更新提供的颜色键，保留未提供的默认值
            for key, value in colors.items():
                if key in self.colors:
                    self.colors[key] = value
                else:
                    warnings.warn(f"Unknown color key: {key}. Using default value.")
        
        # 存储位置信息
        self.aa_positions = []
        self.row_y_positions = {}
        self.all_contiguous_regions = []
        self.beta_contiguous_regions = []
        self.sequence = None
        self.secondary_structure = None
        self.fig = None
        self.ax = None
        self.highlight_regions = []  # 存储高亮区域信息
        self.disulfide_bonds = []  # 存储二硫键连接 [(cys1_index, cys2_index), ...]
    
    def set_colors(self, colors: Dict[str, str]):
        """设置自定义颜色方案"""
        self.colors = {**self.DEFAULT_COLORS, **colors}
    
    
    def set_disulfide_bonds(self, bonds: List[Tuple[int, int]]):
        """
        设置二硫键连接
        
        Args:
            bonds: 二硫键连接列表，每个元素是包含两个半胱氨酸残基索引的元组 (1-based索引)
        """
        self.disulfide_bonds = bonds
    
    def add_highlight_region(self, start: int, end: int, 
                            aa_background: Optional[str] = None, 
                            ss_background: Optional[str] = None,
                            title: Optional[str] = None,
                            title_color: Optional[str] = None,
                            title_fontsize: int = 12):
        """
        添加高亮区域
        
        Args:
            start: 起始氨基酸位置 (1-based)
            end: 结束氨基酸位置 (1-based)
            aa_background: 氨基酸背景颜色
            ss_background: 二级结构背景颜色
            title: 区域标题
            title_color: 标题颜色
            title_fontsize: 标题字体大小
        """
        # 转换为0-based索引
        start_idx = max(0, start - 1)
        end_idx = min(len(self.sequence) - 1, end - 1) if self.sequence else end - 1
        
        self.highlight_regions.append({
            'start': start_idx,
            'end': end_idx,
            'aa_background': aa_background or self.colors.get('highlight_aa', '#FFE0B2'),
            'ss_background': ss_background or self.colors.get('highlight_ss', '#E1F5FE'),
            'title': title,
            'title_color': title_color or self.colors.get('highlight_text', '#D32F2F'),
            'title_fontsize': title_fontsize
        })
    
    def clear_highlight_regions(self):
        """清除所有高亮区域"""
        self.highlight_regions = []
    
    def create_figure(self, 
                 sequence: str,
                 secondary_structure: str,
                 amino_acids_per_line: Optional[int] = None,
                 show_sequence: bool = True,
                 show_residue_numbers: bool = False,
                 show_legend: bool = True,
                 legend_position: str = 'right',
                 title: Optional[str] = "Protein Secondary Structure",
                 dpi: int = 100,
                 figsize: Optional[Tuple[float, float]] = None,
                 adjust_layout: bool = True,
                 show_disulfide_bonds: bool = True,
                 disulfide_bonds: Optional[List[Tuple[int, int]]] = None,
                 ss_range: Optional[Tuple[int, int]] = None,  # 新增参数：二级结构在序列上的位置范围
                 fill_gaps: bool = False,  # 新增参数：控制未覆盖区域是否显示为coil
                 **kwargs) -> plt.Figure:
        """
        Create professional scientific-style secondary structure visualization
        
        Args:
            sequence: Amino acid sequence
            secondary_structure: Secondary structure prediction/assignment
            amino_acids_per_line: Number of amino acids per line (uses instance value if None)
            show_sequence: Whether to show amino acid sequence
            show_residue_numbers: 是否显示每个氨基酸的序号
            show_legend: 是否显示图例
            legend_position: 图例位置 ('bottom' 或 'right')
            title: Figure title
            dpi: Resolution for output
            figsize: 自定义图形尺寸 (宽度, 高度)
            adjust_layout: 是否自动调整布局
            show_disulfide_bonds: 是否显示二硫键连接
            disulfide_bonds: 二硫键连接列表 [(cys1_index, cys2_index), ...] (1-based索引)
            ss_range: 二级结构在序列上的位置范围 (start, end)，1-based索引
            fill_gaps: 是否填充未覆盖区域为coil
            **kwargs: Additional visualization parameters
            
        Returns:
            matplotlib Figure object
        """
        # 处理二级结构范围
        if ss_range is None:
            # 默认情况下，序列和二级结构长度必须相同
            if len(sequence) != len(secondary_structure):
                raise ValueError("Sequence and secondary structure must have same length when ss_range is not provided")
            ss_start = 0
            ss_end = len(secondary_structure) - 1
        else:
            ss_start, ss_end = ss_range
            # 转换为0-based索引
            ss_start_idx = max(0, ss_start - 1)
            ss_end_idx = min(len(secondary_structure) - 1, ss_end - 1)
            
            # 检查范围有效性
            if ss_start_idx > ss_end_idx:
                raise ValueError("Invalid ss_range: start must be less than or equal to end")
            if ss_end_idx - ss_start_idx + 1 > len(sequence):
                raise ValueError("ss_range length exceeds sequence length")
            
            # 调整二级结构字符串以匹配序列长度
            if fill_gaps:
                adjusted_ss = ['C'] * len(sequence)
            else:
                adjusted_ss = [''] * len(sequence)
            for i in range(ss_start_idx, ss_end_idx + 1):
                seq_idx = i - ss_start_idx
                if seq_idx < len(sequence):
                    adjusted_ss[seq_idx] = secondary_structure[i]
            secondary_structure = ''.join(adjusted_ss)
            ss_start = 0
            ss_end = len(secondary_structure) - 1
        
        self.sequence = sequence
        self.secondary_structure = secondary_structure
        
        # 使用实例row_length或覆盖参数
        if amino_acids_per_line is not None:
            self.row_length = amino_acids_per_line
        
        self._precompute_positions(len(sequence))
        self._identify_contiguous_regions(secondary_structure)
        
        total_rows = len(self.row_y_positions)
        
        # 计算图形大小 - 考虑比例因子
        width_ratio = self.row_length / 50
        fig_width = self.base_fig_width * width_ratio * 1.15
        fig_width = min(fig_width, 120) * self.fig_width_scale
        fig_height = total_rows * self.base_fig_height_per_row * self.fig_height_scale
        
        # 如果提供了自定义figsize，则使用
        if figsize is not None:
            fig_width, fig_height = figsize
        
        # 应用宽高比
        if self.aspect_ratio is not None:
            if isinstance(self.aspect_ratio, (int, float)):
                fig_height = fig_width / self.aspect_ratio
            elif isinstance(self.aspect_ratio, str) and ':' in self.aspect_ratio:
                width_ratio, height_ratio = map(float, self.aspect_ratio.split(':'))
                fig_height = fig_width * height_ratio / width_ratio
        
        # 创建图形
        self.fig, self.ax = plt.subplots(figsize=(fig_width, fig_height), dpi=dpi)
        self.ax.set_facecolor(self.colors['background'])
        
        # 绘制高亮区域背景
        self._draw_highlight_regions()
        
        # 绘制所有二级结构
        for (start, end, struct_type) in self.all_contiguous_regions:
            start_row = start // self.row_length
            end_row = end // self.row_length
            
            if start_row == end_row:
                if struct_type == 'H':
                    self._draw_helix_segment(self.ax, start, end, 'H')
                elif struct_type == 'G':
                    self._draw_helix_segment(self.ax, start, end, 'G')
                elif struct_type == 'E':
                    self._draw_beta_strand(self.ax, start, end)
                elif struct_type == 'S':
                    self._draw_bend(self.ax, start, end)
                elif struct_type == 'T':
                    self._draw_turn(self.ax, start, end)
                else:
                    self._draw_coil(self.ax, start, end)
            else:
                current = start
                for row in range(start_row, end_row + 1):
                    row_end = min((row + 1) * self.row_length - 1, end)
                    if current > row_end:
                        continue
                    
                    if struct_type == 'H':
                        self._draw_helix_segment(self.ax, current, row_end, 'H')
                    elif struct_type == 'G':
                        self._draw_helix_segment(self.ax, current, row_end, 'G')
                    elif struct_type == 'I':
                        self._draw_helix_segment(self.ax, current, row_end, 'I')
                    elif struct_type == 'E':
                        self._draw_beta_strand(self.ax, current, row_end)
                    elif struct_type == 'S':
                        self._draw_bend(self.ax, current, row_end)
                    elif struct_type == 'T':
                        self._draw_turn(self.ax, current, row_end)
                    else:
                        self._draw_coil(self.ax, current, row_end)
                    
                    current = row_end + 1
        
        # 绘制氨基酸和残基编号
        if show_sequence:
            self._draw_amino_acids(self.ax, show_residue_numbers=show_residue_numbers)
        
        # 绘制高亮区域标题
        self._draw_highlight_titles()
        
        # 计算边界
        min_x = min(x for x, _ in self.aa_positions) - 2*self.aa_spacing - self.max_element_extension['left']
        max_x = max(x for x, _ in self.aa_positions) + 2*self.aa_spacing + self.max_element_extension['right']
        min_y = min(self.row_y_positions.values()) - self.max_element_extension['bottom']
        max_y = max(self.row_y_positions.values()) + self.max_element_extension['top']
        
        # 设置标题和边界
        self.ax.set_title(title, fontsize=14 + (self.aa_font_size - 10)*0.5, 
                        fontweight='bold', pad=15, color='#5D534A',
                        fontfamily=self.aa_font_family)
        self.ax.set_xlim(min_x, max_x)
        self.ax.set_ylim(min_y, max_y)
        self.ax.axis('off')
        
        # 绘制图例
        if show_legend:
            self._draw_legend(legend_position)  # 传递图例位置参数
        
        # 调整布局
        if adjust_layout:
            self._adjust_layout()

        if show_disulfide_bonds and self.disulfide_bonds:
            self._draw_disulfide_bonds(self.ax)
        
        # 设置二硫键
        if disulfide_bonds is not None:
            self.disulfide_bonds = disulfide_bonds

        # 绘制二硫键连接
        if show_disulfide_bonds and self.disulfide_bonds:
            self._draw_disulfide_bonds(self.ax)
        
        return self.fig

    def create_alignment_figure(self,
                                alignment: Mapping[str, str],
                                structure_sequence: str,
                                secondary_structure: str,
                                amino_acids_per_line: Optional[int] = None,
                                show_legend: bool = False,
                                legend_position: str = 'right',
                                title: Optional[str] = None,
                                dpi: int = 300,
                                figsize: Optional[Tuple[float, float]] = None,
                                structure_label: Optional[str] = None,
                                structure_mode: str = 'major',
                                sequence_colors: Optional[Mapping[str, str]] = None,
                                disulfide_bonds: Optional[Sequence[Tuple[int, int]]] = None,
                                show_disulfide_bonds: bool = True,
                                epitope_annotations: Optional[Mapping[str, Sequence[int]]] = None,
                                epitope_colors: Optional[Mapping[str, str]] = None,
                                show_epitope_legend: bool = True,
                                epitope_marker: str = '^',
                                epitope_marker_size: float = 13.0,
                                epitope_legend_marker_size: float = 4.4,
                                epitope_legend_position: str = 'lower right',
                                epitope_legend_columns: Optional[int] = None,
                                sequence_font_size: Optional[float] = None,
                                label_font_size: Optional[float] = None,
                                sequence_pitch: float = 0.18,
                                structure_offset: float = 0.15,
                                annotation_offset: float = 0.12,
                                annotation_track_spacing: float = 0.105,
                                structure_gap: float = 0.08,
                                show_residue_numbers: bool = True,
                                show_conservation: bool = True,
                                conservation_threshold: float = 0.75,
                                conservation_full_color: str = '#EAF0EE',
                                conservation_partial_color: str = '#F3F5F4',
                                disulfide_color: str = '#2E9663',
                                disulfide_text_color: str = '#237A51',
                                disulfide_line_width: float = 0.65,
                                disulfide_line_length: float = 0.28,
                                disulfide_font_size: float = 5.2,
                                disulfide_line_offset: float = 0.05,
                                disulfide_label_offset: float = -0.05,
                                legend_columns: int = 4,
                                legend_font_size: Optional[float] = None) -> plt.Figure:
        """Draw a compact, publication-style alignment with one structure track.

        The alignment uses a dedicated coordinate system: one x unit is one
        alignment column and the structure glyphs occupy a short annotation band
        above the sequence rows.  This deliberately differs from the spacious
        single-sequence canvas used by :meth:`create_figure`. ``structure_label``
        defaults to the selected sequence name and may be replaced with any text.
        ``structure_mode='major'`` (the default) renders helices, beta strands,
        and turns; ``'all'`` also renders bends and coil segments. Disulfide and
        epitope positions use 1-based, ungapped coordinates from the selected
        structure sequence. The alignment-specific style controls below keep
        the current publication style as defaults while allowing fine tuning of
        typography, spacing, legends, conservation bands, and annotations.

        Common customization controls:
            ``sequence_font_size``, ``label_font_size``, ``sequence_pitch`` and
            ``structure_offset`` control typography and row spacing;
            ``show_residue_numbers``, ``show_conservation`` and
            ``conservation_threshold`` control alignment scaffolding;
            ``legend_position``, ``legend_columns`` and ``legend_font_size``
            control the structure legend; ``show_epitope_legend``,
            ``epitope_marker``, ``epitope_marker_size``,
            ``epitope_legend_position`` and ``epitope_legend_columns`` control
            epitope markers and their legend; ``disulfide_color``,
            ``disulfide_text_color``, ``disulfide_line_width``,
            ``disulfide_line_length``, ``disulfide_font_size``,
            ``disulfide_line_offset`` and ``disulfide_label_offset`` control
            the paired disulfide annotations. Structure colors remain
            configurable through the visualizer ``colors`` argument.
        """
        if not alignment:
            raise ValueError("Alignment cannot be empty")
        if not isinstance(alignment, Mapping):
            raise TypeError("alignment must be a mapping of names to sequences")
        if not isinstance(secondary_structure, str):
            raise TypeError("secondary_structure must be a string")
        if not isinstance(structure_sequence, str):
            raise TypeError("structure_sequence must be a sequence name")
        if structure_label is not None and not isinstance(structure_label, str):
            raise TypeError("structure_label must be a string or None")
        if structure_mode not in {'major', 'all'}:
            raise ValueError("structure_mode must be 'major' or 'all'")
        if sequence_colors is not None and (
                not isinstance(sequence_colors, Mapping)
                or any(not isinstance(name, str) or not isinstance(color, str)
                       for name, color in sequence_colors.items())):
            raise TypeError("sequence_colors must map sequence names to color strings")
        if not isinstance(show_disulfide_bonds, bool):
            raise TypeError("show_disulfide_bonds must be a boolean")
        if not isinstance(show_epitope_legend, bool) or not isinstance(show_residue_numbers, bool):
            raise TypeError("show_epitope_legend and show_residue_numbers must be booleans")
        if not isinstance(show_conservation, bool):
            raise TypeError("show_conservation must be a boolean")
        if epitope_annotations is not None and not isinstance(epitope_annotations, Mapping):
            raise TypeError("epitope_annotations must map labels to residue positions")
        if epitope_colors is not None and (
                not isinstance(epitope_colors, Mapping)
                or any(not isinstance(name, str) or not isinstance(color, str)
                       for name, color in epitope_colors.items())):
            raise TypeError("epitope_colors must map annotation labels to color strings")
        positive_values = {
            'epitope_marker_size': epitope_marker_size,
            'epitope_legend_marker_size': epitope_legend_marker_size,
            'sequence_pitch': sequence_pitch,
            'annotation_track_spacing': annotation_track_spacing,
            'disulfide_line_width': disulfide_line_width,
            'disulfide_line_length': disulfide_line_length,
            'disulfide_font_size': disulfide_font_size,
        }
        if sequence_font_size is not None:
            positive_values['sequence_font_size'] = sequence_font_size
        if label_font_size is not None:
            positive_values['label_font_size'] = label_font_size
        if legend_font_size is not None:
            positive_values['legend_font_size'] = legend_font_size
        if any(not isinstance(value, (int, float)) or value <= 0
               for value in positive_values.values()):
            raise ValueError("font sizes, spacing, and line dimensions must be positive")
        if not isinstance(structure_gap, (int, float)) or structure_gap < 0:
            raise ValueError("structure_gap must be a non-negative number")
        offsets = (structure_offset, annotation_offset,
                   disulfide_line_offset, disulfide_label_offset)
        if any(not isinstance(value, (int, float)) for value in offsets):
            raise TypeError("structure and annotation offsets must be numeric")
        if not isinstance(conservation_threshold, (int, float)) \
                or not 0 < conservation_threshold <= 1:
            raise ValueError("conservation_threshold must be in the interval (0, 1]")
        if not isinstance(legend_columns, int) or legend_columns <= 0:
            raise ValueError("legend_columns must be a positive integer")
        if epitope_legend_columns is not None \
                and (not isinstance(epitope_legend_columns, int) or epitope_legend_columns <= 0):
            raise ValueError("epitope_legend_columns must be a positive integer or None")
        if not isinstance(epitope_marker, str) or not epitope_marker:
            raise TypeError("epitope_marker must be a non-empty Matplotlib marker string")
        if not isinstance(legend_position, str) or not isinstance(epitope_legend_position, str):
            raise TypeError("legend positions must be strings")
        if annotation_offset < 0 or annotation_track_spacing < \
                disulfide_line_offset - disulfide_label_offset:
            raise ValueError(
                "annotation_offset must be non-negative and annotation_track_spacing "
                "must leave room between disulfide lines and labels"
            )

        items = list(alignment.items())
        if any(not isinstance(name, str) or not name or not isinstance(sequence, str)
               for name, sequence in items):
            raise TypeError("Alignment must map non-empty sequence names to strings")

        alignment_length = len(items[0][1])
        if alignment_length == 0 or any(len(sequence) != alignment_length
                                        for _, sequence in items):
            raise ValueError("All aligned sequences must have the same non-zero length")
        if structure_sequence not in alignment:
            raise ValueError(f"Unknown structure sequence: {structure_sequence}")

        per_line = self.row_length if amino_acids_per_line is None else amino_acids_per_line
        if not isinstance(per_line, int) or per_line <= 0:
            raise ValueError("amino_acids_per_line must be a positive integer")

        selected = alignment[structure_sequence]
        ungapped_length = sum(char not in '-.' for char in selected)
        raw_bonds = self.disulfide_bonds if disulfide_bonds is None else disulfide_bonds
        try:
            bonds = [tuple(bond) for bond in raw_bonds]
        except TypeError as exc:
            raise TypeError("disulfide_bonds must contain residue-position pairs") from exc
        if any(len(bond) != 2 or any(not isinstance(position, int) for position in bond)
               for bond in bonds):
            raise TypeError("disulfide_bonds must contain integer residue-position pairs")

        epitope_items = []
        for name, positions in (epitope_annotations or {}).items():
            if not isinstance(name, str) or not name or isinstance(positions, (str, bytes)):
                raise TypeError("epitope_annotations must map labels to integer positions")
            try:
                positions = tuple(positions)
            except TypeError as exc:
                raise TypeError(
                    "epitope_annotations must map labels to integer positions"
                ) from exc
            if any(not isinstance(position, int) for position in positions):
                raise TypeError("epitope annotations must contain integer positions")
            epitope_items.append((name, positions))

        annotated_positions = [position for bond in bonds for position in bond]
        annotated_positions.extend(
            position for _, positions in epitope_items for position in positions
        )
        if any(position < 1 or position > ungapped_length
               for position in annotated_positions):
            raise ValueError(
                "Disulfide and epitope positions must fall within the selected sequence"
            )
        residue_columns = {}
        residue_position = 0
        for column, residue in enumerate(selected):
            if residue not in '-.':
                residue_position += 1
                residue_columns[residue_position] = column

        secondary_structure = secondary_structure.upper()
        if len(secondary_structure) == ungapped_length:
            structure_iter = iter(secondary_structure)
            aligned_structure = ''.join(
                ' ' if char in '-.' else next(structure_iter) for char in selected
            )
        elif len(secondary_structure) == alignment_length:
            aligned_structure = ''.join(
                ' ' if char in '-.' else ss
                for char, ss in zip(selected, secondary_structure)
            )
        else:
            raise ValueError(
                "Secondary structure length must match the selected sequence "
                "with or without alignment gaps"
            )

        aligned_structure = aligned_structure.replace('-', ' ').replace('.', ' ')
        invalid_codes = set(aligned_structure) - set('HGIESTC ')
        if invalid_codes:
            raise ValueError(f"Unknown secondary structure codes: {sorted(invalid_codes)}")

        self.sequence = selected
        self.secondary_structure = aligned_structure
        self.row_length = per_line
        displayed_structure_label = structure_sequence \
            if structure_label is None else structure_label
        block_count = (alignment_length + per_line - 1) // per_line
        display_columns = min(per_line, alignment_length)
        sequence_font_size = (
            min(8.0, max(6.5, self.aa_font_size * 0.48))
            if sequence_font_size is None else float(sequence_font_size)
        )
        label_font_size = (
            max(6.5, sequence_font_size * 0.94)
            if label_font_size is None else float(label_font_size)
        )
        has_disulfide_track = show_disulfide_bonds and bool(bonds)
        annotation_pitch = float(annotation_track_spacing)
        annotation_track_count = len(epitope_items) + int(has_disulfide_track)
        annotation_depth = (
            annotation_offset + (annotation_track_count - 1) * annotation_pitch
            if annotation_track_count else 0.0
        )
        block_height = len(items) * sequence_pitch + 0.55 + annotation_depth
        block_top = {block: -block * block_height for block in range(block_count)}
        self.row_y_positions = {
            block: block_top[block] - structure_offset for block in range(block_count)
        }
        self.aa_positions = [
            (index % per_line + 0.5, index // per_line)
            for index in range(alignment_length)
        ]
        self._identify_contiguous_regions(aligned_structure)

        longest_name = max(len(displayed_structure_label),
                           max(len(name) for name, _ in items))
        label_inches = longest_name * label_font_size * 0.60 / 72
        width_inches = max(
            3.8,
            min(12.5, 0.45 + label_inches + display_columns * 0.105)
        )
        top_space = 0.30 if title or show_legend else 0.12
        height_inches = max(
            1.50,
            top_space + block_count * (
                0.30 + 0.16 * len(items) + 0.75 * annotation_depth
            ) + (0.16 if epitope_items else 0)
        )
        fig_width = width_inches * self.fig_width_scale
        fig_height = height_inches * self.fig_height_scale
        if figsize is None:
            figsize = (fig_width, fig_height)

        self.fig, self.ax = plt.subplots(figsize=figsize, dpi=dpi)
        # The default single-sequence canvas is warm, while alignments are
        # rendered on white to match current structural-biology figures.
        background = '#FFFFFF' if self.colors['background'] == self.DEFAULT_COLORS['background'] \
            else self.colors['background']
        self.fig.set_facecolor(background)
        self.ax.set_facecolor(background)

        structure_colors = {
            'H': '#3E8F5A', 'G': '#B96F91', 'I': '#4F70C8',
            'E': '#E5BF70', 'T': '#9D8BC5', 'S': '#55AFC2',
            'C': '#A7ADB1',
        }
        for code in structure_colors:
            if code in self.colors and self.colors[code] != self.DEFAULT_COLORS[code]:
                structure_colors[code] = self.colors[code]
        visible_codes = set('HGIET') if structure_mode == 'major' else set('HGIESTC')
        sequence_colors = dict(sequence_colors or {})
        default_epitope_colors = {
            'h5C_4 epitope': '#8FB9DB',
            'C05 epitope': '#D39A9B',
            'D07 epitope': '#AFC56F',
        }
        fallback_epitope_colors = ('#5D9DC6', '#C47F82', '#8EAA4E', '#8B78B5')
        resolved_epitope_colors = {}
        for index, (name, _) in enumerate(epitope_items):
            resolved_epitope_colors[name] = (epitope_colors or {}).get(
                name,
                default_epitope_colors.get(
                    name, fallback_epitope_colors[index % len(fallback_epitope_colors)]
                ),
            )
        structure_label_color = self.colors.get('structure_label', '#616971')
        mono_family = 'DejaVu Sans Mono'
        backbone_color = '#B5BABE'

        def draw_structure(ax, start, end, code, y, terminal=True,
                           gap_before=0.0, gap_after=0.0):
            """Draw one structure segment on shared integer column edges."""
            left = float(start) + gap_before
            right = float(end + 1) - gap_after
            if right <= left:
                return
            color = structure_colors[code]
            if code in {'H', 'G', 'I'}:
                # Preserve PyPlogo's original helix language: each residue is
                # a dark/light overlapping parallelogram pair.  The geometry
                # is scaled to one alignment column instead of reusing the
                # large single-sequence dimensions.
                # Keep the original 15:12:12 width/overlap/skew relationship,
                # with a small compactness correction for the tighter MSA row.
                # The reduced height keeps the motif subordinate to the
                # residue glyphs while the brighter facet preserves the
                # interlocking diamond rhythm at journal figure scale.
                light = self._adjust_color(color, 1.42)
                compact = 0.90
                width = compact * 15 / 28
                overlap = compact * 12 / 28
                skew = compact * 12 / 28
                half_height = 0.064
                unit_width = 2 * width - overlap
                x_adjustment = (1.0 - unit_width) / 2
                clip = patches.Rectangle(
                    (left, y - 0.14), right - left, 0.28,
                    transform=ax.transData
                )
                for cell in range(start, end + 1):
                    center = cell + 0.5 + x_adjustment
                    dark_vertices = [
                        (center - width, y - half_height),
                        (center - width - skew, y + half_height),
                        (center - skew, y + half_height),
                        (center, y - half_height),
                    ]
                    light_vertices = [
                        (center - overlap, y - half_height),
                        (center - overlap + skew, y + half_height),
                        (center - overlap + width + skew, y + half_height),
                        (center - overlap + width, y - half_height),
                    ]
                    for vertices, facecolor, zorder in (
                        (light_vertices, light, 3),
                        (dark_vertices, color, 4),
                    ):
                        patch = patches.Polygon(
                            vertices, closed=True, facecolor=facecolor,
                            edgecolor='none', zorder=zorder
                        )
                        patch.set_clip_path(clip)
                        ax.add_patch(patch)
            elif code == 'E':
                width = right - left
                body = 0.034
                if terminal:
                    # Keep a prominent arrowhead inside this beta segment so
                    # neighboring structure elements remain unobscured.
                    head = min(1.25, max(0.35, width * 0.26), width * 0.42)
                    head_half = min(0.115, max(0.085, width * 0.16))
                    points = [
                        (left, y - body), (right - head, y - body),
                        (right - head, y - head_half), (right, y),
                        (right - head, y + head_half), (right - head, y + body),
                        (left, y + body),
                    ]
                else:
                    # A wrapped beta strand continues as a flat ribbon; the
                    # arrowhead appears only at the true C-terminal end.
                    points = [
                        (left, y - body), (right, y - body),
                        (right, y + body), (left, y + body),
                    ]
                ax.add_patch(patches.Polygon(
                    points, closed=True, facecolor=color, edgecolor='#A47B37',
                    linewidth=0.32, joinstyle='miter', zorder=3
                ))
            elif code == 'T':
                turn_label = 'T' * min(3, end - start + 1)
                ax.text(
                    (left + right) / 2, y, turn_label,
                    ha='center', va='center', fontsize=6.1,
                    color=color, fontweight='bold', fontfamily='DejaVu Sans',
                    zorder=4,
                )
            else:
                if code == 'S':
                    bend_path = Path(
                        [(left, y), ((left + right) / 2, y + 0.052), (right, y)],
                        [Path.MOVETO, Path.LINETO, Path.LINETO]
                    )
                    ax.add_patch(patches.PathPatch(
                        bend_path, fill=False, edgecolor=color,
                        linewidth=1.30, capstyle='round', joinstyle='round',
                        zorder=3
                    ))

        label_x = -0.34
        min_x = -max(2.2, longest_name * 0.58) - 0.10
        max_x = display_columns + 0.18
        for block in range(block_count):
            start = block * per_line
            end = min(start + per_line, alignment_length)
            columns = end - start
            number_y = block_top[block]
            ss_y = self.row_y_positions[block]
            sequence_top = number_y - 0.34
            sequence_bottom = sequence_top - (len(items) - 1) * sequence_pitch

            # Faint conservation bands connect related residues vertically
            # without recreating the saturated ESPript color wall.
            for column in range(start, end):
                if not show_conservation:
                    break
                residues = [
                    sequence[column].upper() for _, sequence in items
                    if sequence[column] not in '-.'
                ]
                if len(residues) < 2:
                    continue
                counts = {residue: residues.count(residue) for residue in set(residues)}
                fraction = max(counts.values()) / len(residues)
                if fraction == 1.0 and len(residues) == len(items):
                    band_color = conservation_full_color
                elif fraction >= conservation_threshold:
                    band_color = conservation_partial_color
                else:
                    continue
                local_column = column - start
                self.ax.add_patch(patches.Rectangle(
                    (local_column + 0.04, sequence_bottom - 0.10),
                    0.92, sequence_top - sequence_bottom + 0.20,
                    facecolor=band_color, edgecolor='none', zorder=0
                ))

            # The major mode contains structure glyphs only. Full mode also
            # retains dotted alignment-gap connectors and explicit coil runs.
            if structure_mode == 'all':
                for column, residue in enumerate(selected[start:end]):
                    if residue not in '-.':
                        continue
                    self.ax.plot(
                        [column + 0.12, column + 0.88], [ss_y, ss_y],
                        color=backbone_color, linewidth=0.70,
                        linestyle=(0, (1.0, 1.0)), zorder=2
                    )

            if 'C' in visible_codes:
                coil_start = None
                for local_column in range(columns + 1):
                    is_coil = (
                        local_column < columns
                        and aligned_structure[start + local_column] == 'C'
                    )
                    if is_coil and coil_start is None:
                        coil_start = local_column
                    elif not is_coil and coil_start is not None:
                        self.ax.plot(
                            [coil_start, local_column], [ss_y, ss_y],
                            color=backbone_color, linewidth=0.85,
                            solid_capstyle='butt', zorder=1
                        )
                        coil_start = None

            for region_start, region_end, code in self.all_contiguous_regions:
                if code not in visible_codes or code == 'C':
                    continue
                segment_start = max(region_start, start)
                segment_end = min(region_end, end - 1)
                if segment_start <= segment_end:
                    draw_structure(
                        self.ax, segment_start - start, segment_end - start,
                        code, ss_y,
                        terminal=(code != 'E' or segment_end == region_end),
                        gap_before=(structure_gap
                                    if segment_start == region_start
                                    and segment_start > start
                                    and aligned_structure[region_start - 1] in visible_codes
                                    else 0.0),
                        gap_after=(structure_gap
                                   if segment_end == region_end
                                   and segment_end + 1 < end else 0.0),
                    )

            if show_residue_numbers:
                self.ax.text(
                    0, number_y, str(start + 1), ha='left', va='bottom',
                    fontsize=max(5.8, sequence_font_size - 1.3),
                    color=self.colors['residue_number'], fontweight='normal',
                    fontfamily=mono_family, zorder=5
                )
                self.ax.text(
                    columns, number_y, str(end), ha='right', va='bottom',
                    fontsize=max(5.8, sequence_font_size - 1.3),
                    color=self.colors['residue_number'], fontweight='normal',
                    fontfamily=mono_family, zorder=5
                )
            if displayed_structure_label:
                self.ax.text(
                    label_x, ss_y, displayed_structure_label, ha='right', va='center',
                    fontsize=label_font_size - 0.3, color=structure_label_color,
                    fontfamily='DejaVu Sans', zorder=5
                )

            for row, (name, sequence) in enumerate(items):
                text_y = sequence_top - row * sequence_pitch
                sequence_color = sequence_colors.get(name, self.colors['sequence'])
                self.ax.text(
                    label_x, text_y, name, ha='right', va='center',
                    fontsize=label_font_size, color=sequence_color,
                    fontweight='bold' if name == structure_sequence else 'normal',
                    fontfamily='DejaVu Sans', zorder=5
                )
                for column, residue in enumerate(sequence[start:end]):
                    residue = '-' if residue == '.' else residue
                    self.ax.text(
                        column + 0.5, text_y, residue, ha='center', va='center',
                        fontsize=sequence_font_size,
                        color=self.colors['residue_number'] if residue == '-' else sequence_color,
                        fontweight='bold' if name == structure_sequence else 'normal',
                        fontfamily=mono_family, zorder=5
                    )

            annotation_y = sequence_bottom - annotation_offset
            if has_disulfide_track:
                for bond_id, bond in enumerate(bonds, 1):
                    for position in bond:
                        column = residue_columns[position]
                        if start <= column < end:
                            x = column - start + 0.5
                            marker_line, = self.ax.plot(
                                [x - disulfide_line_length / 2,
                                 x + disulfide_line_length / 2],
                                [annotation_y + disulfide_line_offset,
                                 annotation_y + disulfide_line_offset],
                                color=disulfide_color, linewidth=disulfide_line_width,
                                solid_capstyle='round', zorder=5,
                            )
                            marker_line.set_gid(f'disulfide-{bond_id}')
                            marker = self.ax.text(
                                x, annotation_y + disulfide_label_offset, str(bond_id),
                                ha='center', va='center', fontsize=disulfide_font_size,
                                color=disulfide_text_color, fontweight='bold',
                                fontfamily='DejaVu Sans', zorder=6,
                            )
                            marker.set_gid(f'disulfide-{bond_id}')
                annotation_y -= annotation_pitch

            for name, positions in epitope_items:
                x_positions = [
                    residue_columns[position] - start + 0.5
                    for position in positions
                    if start <= residue_columns[position] < end
                ]
                if x_positions:
                    markers = self.ax.scatter(
                        x_positions, [annotation_y] * len(x_positions),
                        marker=epitope_marker, s=epitope_marker_size,
                        color=resolved_epitope_colors[name],
                        edgecolors='none', linewidths=0, zorder=5,
                    )
                    markers.set_gid(f'epitope-{name}')
                annotation_y -= annotation_pitch

        min_y = -(block_count - 1) * block_height \
            - 0.34 - (len(items) - 1) * sequence_pitch - 0.13 - annotation_depth
        max_y = 0.07
        self.ax.set_xlim(min_x, max_x)
        self.ax.set_ylim(min_y, max_y)
        self.ax.axis('off')
        if title:
            self.ax.set_title(
                title, loc='left', fontsize=9.0, fontweight='bold',
                color='#27313A', fontfamily='DejaVu Sans', pad=7
            )

        if show_legend:
            from matplotlib.lines import Line2D
            handles = []
            if 'H' in visible_codes:
                handles.append(patches.Patch(
                    facecolor=structure_colors['H'], edgecolor='none', label='alpha-helix'))
            if 'G' in visible_codes:
                handles.append(patches.Patch(
                    facecolor=structure_colors['G'], edgecolor='none', label='3_10-helix'))
            if 'I' in visible_codes:
                handles.append(patches.Patch(
                    facecolor=structure_colors['I'], edgecolor='none', label='pi-helix'))
            if 'E' in visible_codes:
                handles.append(patches.Patch(
                    facecolor=structure_colors['E'], edgecolor='#A47B37', label='beta-strand'))
            if 'T' in visible_codes:
                handles.append(Line2D(
                    [0], [0], color='none', marker='$TT$', markersize=7,
                    markeredgecolor=structure_colors['T'], label='turn'
                ))
            if 'S' in visible_codes:
                handles.append(Line2D([0], [0], color=structure_colors['S'], lw=1.5, label='bend'))
            if 'C' in visible_codes:
                handles.append(Line2D([0], [0], color=structure_colors['C'], lw=1.2, label='coil'))
            legend_locations = {
                'right': ('upper right', (0.995, 0.995)),
                'left': ('upper left', (0.005, 0.995)),
                'top': ('upper center', (0.50, 0.995)),
                'bottom': ('lower center', (0.50, 0.005)),
            }
            legend_loc, legend_anchor = legend_locations.get(
                legend_position, (legend_position, None)
            )
            self.fig.legend(
                handles=handles, ncol=legend_columns, loc=legend_loc,
                bbox_to_anchor=legend_anchor, frameon=False,
                fontsize=(max(6.0, sequence_font_size - 1.1)
                          if legend_font_size is None else legend_font_size),
                handlelength=1.15, handletextpad=0.30, columnspacing=0.65,
            )
        if epitope_items and show_epitope_legend:
            from matplotlib.lines import Line2D
            epitope_handles = [
                Line2D(
                    [0], [0], linestyle='none', marker=epitope_marker,
                    markersize=epitope_legend_marker_size,
                    markerfacecolor=resolved_epitope_colors[name],
                    markeredgewidth=0, label=name,
                )
                for name, _ in epitope_items
            ]
            self.fig.legend(
                handles=epitope_handles,
                ncol=(min(3, len(epitope_handles)) if epitope_legend_columns is None
                      else epitope_legend_columns),
                loc=epitope_legend_position, frameon=False,
                fontsize=(max(6.0, sequence_font_size - 1.0)
                          if legend_font_size is None else legend_font_size),
                handlelength=0.8, handletextpad=0.25, columnspacing=0.75,
            )
        top = 0.82 if title or show_legend else 0.96
        bottom = 0.12 if epitope_items and show_epitope_legend else 0.06
        self.fig.subplots_adjust(left=0.02, right=0.99, top=top, bottom=bottom)

        return self.fig
    
    def _draw_disulfide_bonds(self, ax):
        """绘制二硫键连接标记"""
        # 为每个二硫键分配唯一ID
        bond_id = 1
        self.disulfide_markers = {}
        
        for bond in self.disulfide_bonds:
            cys1, cys2 = bond
            
            # 转换为0-based索引
            cys1_idx = cys1 - 1
            cys2_idx = cys2 - 1
            
            # 确保索引在序列范围内
            if cys1_idx < 0 or cys1_idx >= len(self.sequence) or \
            cys2_idx < 0 or cys2_idx >= len(self.sequence):
                continue
            
            # 获取位置
            if cys1_idx < len(self.aa_positions) and cys2_idx < len(self.aa_positions):
                x1, row1 = self.aa_positions[cys1_idx]
                x2, row2 = self.aa_positions[cys2_idx]
                
                y1 = self.row_y_positions[row1]
                y2 = self.row_y_positions[row2]
                
                # 计算标记位置 (在二级结构下方)
                marker_y1 = y1 - self.helix_unit_height + self.disulfide_height_offset
                marker_y2 = y2 - self.helix_unit_height + self.disulfide_height_offset
                
                # 绘制标记
                self._draw_disulfide_marker(ax, x1, marker_y1, bond_id)
                self._draw_disulfide_marker(ax, x2, marker_y2, bond_id)
                
                # 存储标记信息
                self.disulfide_markers[bond_id] = [cys1_idx, cys2_idx]
                bond_id += 1

    def _draw_disulfide_marker(self, ax, x, y, bond_id):
        """绘制单个二硫键标记"""
        text_color = self.colors['disulfide_text']
        
        # 绘制数字
        ax.text(
            x, y, str(bond_id),
            ha='center', va='center',
            fontsize=self.disulfide_font_size-3,
            color=text_color,
            fontweight='bold',
            fontfamily=self.aa_font_family,
            zorder=11
        )
    
    def save_figure(self, file_path: str, **kwargs):
        """Save current figure to file"""
        if self.fig is None:
            raise ValueError("No figure to save. Call create_figure first.")
        
        self.fig.savefig(file_path, 
                        bbox_inches='tight',
                        facecolor=self.colors['background'],
                        **kwargs)
    
    def _precompute_positions(self, sequence_length: int):
        """Precompute amino acid positions and row Y coordinates"""
        self.aa_positions = []
        self.row_y_positions = {}
        
        total_rows = (sequence_length + self.row_length - 1) // self.row_length
        total_height = total_rows * self.row_height
        start_y = 100 + (total_height // 2)
        
        for row in range(total_rows):
            self.row_y_positions[row] = start_y - (row * self.row_height)
        
        for i in range(sequence_length):
            row = i // self.row_length
            pos_in_row = i % self.row_length
            
            x_center = 100 + pos_in_row * self.aa_spacing + self.aa_spacing / 2
            self.aa_positions.append((x_center, row))
    
    def _identify_contiguous_regions(self, secondary_structure: str):
        """Identify contiguous secondary structure regions"""
        self.all_contiguous_regions = []
        if not secondary_structure:
            return
            
        current_type = secondary_structure[0]
        start = 0
        
        for i in range(1, len(secondary_structure)):
            if secondary_structure[i] != current_type:
                self.all_contiguous_regions.append((start, i-1, current_type))
                current_type = secondary_structure[i]
                start = i
        
        self.all_contiguous_regions.append((start, len(secondary_structure)-1, current_type))
        
        # Extract beta strand regions
        self.beta_contiguous_regions = []
        for (start, end, struct_type) in self.all_contiguous_regions:
            if struct_type == 'E':
                start_row = start // self.row_length
                end_row = end // self.row_length
                self.beta_contiguous_regions.append((start, end, start_row, end_row))
    
    def _draw_highlight_regions(self):
        """绘制高亮区域背景"""
        for region in self.highlight_regions:
            start = region['start']
            end = region['end']
            aa_bg = region['aa_background']
            ss_bg = region['ss_background']
            
            # 找到区域所在的所有行
            start_row = start // self.row_length
            end_row = end // self.row_length
            
            for row in range(start_row, end_row + 1):
                # 计算该行的起始和结束索引
                row_start = max(start, row * self.row_length)
                row_end = min(end, (row + 1) * self.row_length - 1)
                
                if row_start > row_end:
                    continue
                    
                # 获取该行的Y坐标
                y_center = self.row_y_positions[row]
                
                # 计算氨基酸序列的位置
                aa_y = y_center + self.helix_unit_height + 5 + (self.aa_font_size - 10) * 0.3
                aa_height = self.aa_font_size * 1.2
                
                # 计算二级结构的位置
                ss_y = y_center
                ss_height = self.helix_unit_height * 2
                
                # 计算该行区域的X范围
                min_x = min(self.aa_positions[i][0] for i in range(row_start, row_end + 1)) - self.aa_spacing / 2
                max_x = max(self.aa_positions[i][0] for i in range(row_start, row_end + 1)) + self.aa_spacing / 2
                width = max_x - min_x
                
                # 绘制氨基酸背景高亮
                if aa_bg:
                    rect = patches.Rectangle(
                        (min_x, aa_y - aa_height / 2),
                        width,
                        aa_height,
                        facecolor=aa_bg,
                        alpha=0.3,
                        edgecolor='none',
                        zorder=0
                    )
                    self.ax.add_patch(rect)
                
                # 绘制二级结构背景高亮
                if ss_bg:
                    rect = patches.Rectangle(
                        (min_x, ss_y - ss_height / 2),
                        width,
                        ss_height,
                        facecolor=ss_bg,
                        alpha=0.2,
                        edgecolor='none',
                        zorder=0
                    )
                    self.ax.add_patch(rect)
    
    def _draw_highlight_titles(self):
        """绘制高亮区域标题 - 修复跨行位置问题"""
        for region in self.highlight_regions:
            if not region.get('title'):
                continue
                
            start = region['start']
            end = region['end']
            title = region['title']
            color = region['title_color']
            fontsize = region['title_fontsize']
            
            # 找到区域所在的所有行
            start_row = start // self.row_length
            end_row = end // self.row_length
            
            # 计算每行的长度
            row_lengths = {}
            for row in range(start_row, end_row + 1):
                row_start = max(start, row * self.row_length)
                row_end = min(end, (row + 1) * self.row_length - 1)
                row_length = row_end - row_start + 1
                row_lengths[row] = row_length
            
            # 找到最长的行
            max_length = 0
            max_row = start_row
            for row, length in row_lengths.items():
                if length > max_length:
                    max_length = length
                    max_row = row
            
            # 计算该行的起始和结束索引
            row_start = max(start, max_row * self.row_length)
            row_end = min(end, (max_row + 1) * self.row_length - 1)
            
            # 计算该行区域的中心位置
            start_x = self.aa_positions[row_start][0]
            end_x = self.aa_positions[row_end][0]
            center_x = (start_x + end_x) / 2
            
            # 计算垂直位置
            y_pos = self.row_y_positions[max_row] + self.helix_unit_height * 2.5
            
            self.ax.text(
                center_x, y_pos, title,
                fontsize=fontsize,
                color=color,
                fontweight='bold',
                ha='center',
                va='bottom',
                fontfamily=self.aa_font_family,
                zorder=20
            )
    
    def _draw_amino_acids(self, ax, show_residue_numbers=False):
        """绘制氨基酸和残基编号 - 去除白色背景框"""
        total_aa = len(self.sequence)
        total_rows = len(self.row_y_positions)
        
        # 计算残基编号的垂直位置偏移量
        # 使用新参数residue_number_position控制位置 (0.0=最低, 1.0=最高)
        residue_y_offset = self.helix_unit_height + 5 + (self.aa_font_size - 10) * 0.3
        residue_y_offset += self.aa_font_size * self.residue_number_position * 1.5
        
        for row in range(total_rows):
            # 计算该行的起始和结束索引
            row_start_idx = row * self.row_length
            row_end_idx = min((row + 1) * self.row_length - 1, total_aa - 1)
            row_y_center = self.row_y_positions[row]
            
            # 绘制起始残基编号
            if not show_residue_numbers:
                first_aa_x = self.aa_positions[row_start_idx][0]
                ax.text(
                    first_aa_x - self.aa_spacing / 2 - 10,
                    row_y_center + residue_y_offset,
                    str(row_start_idx + 1),
                    ha='right', va='bottom',
                    fontsize=self.residue_number_font_size,
                    color=self.colors['residue_number'],
                    fontweight='bold',
                    fontfamily=self.aa_font_family,
                    zorder=10
                )
                
                # 绘制结束残基编号
                last_aa_x = self.aa_positions[row_end_idx][0]
                ax.text(
                    last_aa_x + self.aa_spacing / 2 + 10,
                    row_y_center + residue_y_offset,
                    str(row_end_idx + 1),
                    ha='left', va='bottom',
                    fontsize=self.residue_number_font_size,
                    color=self.colors['residue_number'],
                    fontweight='bold',
                    fontfamily=self.aa_font_family,
                    zorder=10
                )
        
        # 绘制氨基酸序列 - 去除白色背景框
        for i, (x_center, row) in enumerate(self.aa_positions):
            y_center = self.row_y_positions[row]
            text_y = y_center + self.helix_unit_height + 5 + (self.aa_font_size - 10)*0.3
            
            # 绘制氨基酸序号
            if show_residue_numbers:
                residue_y = text_y + self.aa_font_size * self.residue_number_position * 1.5
                ax.text(
                    x_center, residue_y, str(i+1),
                    ha='center', va='bottom',
                    fontsize=self.residue_number_font_size,
                    color=self.colors['residue_number'],
                    fontfamily=self.aa_font_family,
                    zorder=10
                )
            
            # 绘制氨基酸字母 - 去除白色背景框
            ax.text(
                x_center, text_y, self.sequence[i], 
                ha='center', va='center', 
                fontsize=self.aa_font_size,
                color=self.colors['sequence'],
                fontfamily=self.aa_font_family,
                zorder=10
            )
    
    def _draw_legend(self, legend_position: str = 'right'):
        """在指定位置绘制图例"""
        # 图例元素
        elements = [
            ('α-Helix', 'H', self.colors['H'], 'block'),
            ('β-Strand', 'E', self.colors['E'], 'block'),
            ('3₁₀-Helix', 'G', self.colors['G'], 'block'),
            ('π-Helix', 'I', self.colors['I'], 'block'),
            ('Turn', 'T', self.colors['T'], 'symbol'),
            ('Bend', 'S', self.colors['S'], 'symbol'),
            ('Coil', 'C', self.colors['C'], 'symbol')
        ]
        
        if legend_position == 'bottom':
            # 底部图例 - 分两行显示
            self._draw_bottom_legend(elements)
        else:
            # 右侧图例 - 单列显示
            self._draw_right_legend(elements)

    def _draw_bottom_legend(self, elements):
        """在底部绘制图例，一行显示所有元素"""
        # 计算图例位置
        min_x = min(x for x, _ in self.aa_positions) - 2 * self.aa_spacing
        max_x = max(x for x, _ in self.aa_positions) + 2 * self.aa_spacing
        min_y = min(self.row_y_positions.values()) - self.max_element_extension['bottom']
        
        legend_y = min_y + 20
        legend_height = 20
        
        # 计算图例项宽度
        total_width = max_x - min_x
        num_elements = len(elements)
        element_width = total_width / (num_elements + 4)
        
        # 绘制图例 - 所有元素在一行
        for i, (label, ss_type, color, display_type) in enumerate(elements):
            center_x = min_x + element_width * (i + 2.5)
            y_pos = legend_y
            
            if display_type == 'block':
                # 绘制颜色块 (用于螺旋和折叠)
                rect = patches.Rectangle(
                    (center_x - 12, y_pos - 15), 25, 8,
                    facecolor=color,
                    edgecolor=self.colors['edge'],
                    linewidth=0.5,
                    zorder=10
                )
                self.ax.add_patch(rect)
            else:
                # 绘制符号 (用于Turn、Bend和Coil)
                if ss_type == 'T':  # Turn - 弧线
                    x_vals = np.linspace(center_x - 13, center_x + 13, 50)
                    y_vals = y_pos + 3 * np.sin(np.pi * (x_vals - (center_x - 15)) / 30) - 10
                    self.ax.plot(x_vals, y_vals, color=color, linewidth=2.5, zorder=10)
                elif ss_type == 'S':  # Bend - 直线
                    self.ax.plot([center_x - 10, center_x + 10], [y_pos - 10, y_pos - 10], 
                                color=color, linewidth=2.5, zorder=10)
                else:  # Coil - 直线
                    self.ax.plot([center_x - 10, center_x + 10], [y_pos - 10, y_pos - 10], 
                                color=color, linewidth=2.5, zorder=10)
            
            # 绘制标签
            self.ax.text(
                center_x, y_pos - 20, label,
                ha='center', va='top',
                fontsize=self.aa_font_size - 2,
                color=self.colors['legend_text'],
                fontfamily=self.aa_font_family,
                zorder=10
            )

    def _draw_right_legend(self, elements):
        """在右侧绘制图例，单列显示"""
        # 计算图例位置
        max_x = max(x for x, _ in self.aa_positions) + 2 * self.aa_spacing
        min_y = min(self.row_y_positions.values()) - self.max_element_extension['bottom']
        max_y = max(self.row_y_positions.values()) + self.max_element_extension['top']
        
        legend_x = max_x + 40
        legend_y_center = (min_y + max_y) / 2
        
        # 计算图例项高度
        element_height = 20
        total_height = len(elements) * element_height
        start_y = legend_y_center + total_height / 2 - element_height / 2
        
        # 绘制图例
        for i, (label, ss_type, color, display_type) in enumerate(elements):
            y_pos = start_y - i * element_height
            
            if display_type == 'block':
                # 绘制颜色块 (用于螺旋和折叠)
                rect = patches.Rectangle(
                    (legend_x - 12, y_pos - 4), 25, 8,
                    facecolor=color,
                    edgecolor=self.colors['edge'],
                    linewidth=0.5,
                    zorder=10
                )
                self.ax.add_patch(rect)
            else:
                # 绘制符号 (用于Turn、Bend和Coil)
                if ss_type == 'T':  # Turn - 弧线
                    x_vals = np.linspace(legend_x - 13, legend_x + 13, 50)
                    y_vals = y_pos + 3 * np.sin(np.pi * (x_vals - (legend_x - 15)) / 30)
                    self.ax.plot(x_vals, y_vals, color=color, linewidth=2.5, zorder=10)
                elif ss_type == 'S':  # Bend - 直线
                    self.ax.plot([legend_x - 10, legend_x + 10], [y_pos, y_pos], 
                                color=color, linewidth=2.5, zorder=10)
                else:  # Coil - 直线
                    self.ax.plot([legend_x - 10, legend_x + 10], [y_pos, y_pos], 
                                color=color, linewidth=2.5, zorder=10)
            
            # 绘制标签
            self.ax.text(
                legend_x + 25, y_pos, label,
                ha='left', va='center',
                fontsize=self.aa_font_size - 2,
                color=self.colors['legend_text'],
                fontfamily=self.aa_font_family,
                zorder=10
            )
    
    def _draw_helix_symbol(self, x, y, color):
        """绘制螺旋图例符号"""
        width = 20
        height = 10
        skew = 5
        
        # 左平行四边形
        left_vertices = [
            (x - width/2, y - height/2),
            (x - width/2 - skew, y + height/2),
            (x - skew, y + height/2),
            (x, y - height/2)
        ]
        
        # 右平行四边形
        right_vertices = [
            (x - skew, y - height/2),
            (x - skew + skew, y + height/2),
            (x - skew + width + skew, y + height/2),
            (x - skew + width, y - height/2)
        ]
        
        # 绘制
        self._draw_parallelogram(self.ax, left_vertices, color, alpha=0.9)
        self._draw_parallelogram(self.ax, right_vertices, color, alpha=0.7)
    
    def _draw_beta_symbol(self, x, y, color):
        """绘制β折叠图例符号"""
        width = 25
        height = 8
        head_length = 8
        
        # 箭头主体
        body_width = width - head_length
        vertices = [
            (x - width/2, y - height/2),
            (x - width/2 + body_width, y - height/2),
            (x - width/2 + body_width, y - height/4),
            (x - width/2 + width, y),
            (x - width/2 + body_width, y + height/4),
            (x - width/2 + body_width, y + height/2),
            (x - width/2, y + height/2),
            (x - width/2, y - height/2)
        ]
        
        # 绘制
        path = Path(vertices)
        patch = patches.PathPatch(
            path, facecolor=color, alpha=0.9, edgecolor='none', zorder=10
        )
        self.ax.add_patch(patch)
    
    def _draw_ss_symbol(self, x, y, ss_type, color):
        """绘制其他二级结构图例符号"""
        if ss_type == 'T':  # Turn
            # 绘制曲线
            x_vals = np.linspace(x - 15, x + 15, 50)
            y_vals = y + 3 * np.sin(np.pi * (x_vals - (x - 15)) / 30)
            self.ax.plot(x_vals, y_vals, color=color, linewidth=2.5, zorder=10)
        elif ss_type == 'S':  # Bend
            # 绘制直线
            self.ax.plot([x - 15, x + 15], [y, y], color=color, linewidth=2.5, zorder=10)
        else:  # Coil
            # 绘制波浪线
            x_vals = np.linspace(x - 15, x + 15, 50)
            y_vals = y + 2 * np.sin(2 * np.pi * (x_vals - (x - 15)) / 30)
            self.ax.plot(x_vals, y_vals, color=color, linewidth=2.5, zorder=10)
    
    def _draw_parallelogram(self, ax, vertices, color, alpha=1.0, zorder=5):
        """绘制平行四边形"""
        vertices.append(vertices[0])
        codes = [Path.MOVETO] + [Path.LINETO]*(len(vertices)-2) + [Path.CLOSEPOLY]
        path = Path(vertices, codes)
        
        # 填充
        fill_patch = patches.PathPatch(
            path, facecolor=color, alpha=alpha, edgecolor='none', zorder=zorder
        )
        ax.add_patch(fill_patch)
    
    def _draw_helix_segment(self, ax, start, end, helix_type):
        """绘制螺旋段"""
        start_row = start // self.row_length
        end_row = end // self.row_length
        
        # 收集所有氨基酸索引
        aa_indices = []
        if start_row == end_row:
            aa_indices = list(range(start, end + 1))
        else:
            # 第一行
            row_end = ((start_row + 1) * self.row_length) - 1
            aa_indices.extend(range(start, min(row_end, end) + 1))
            
            # 中间行
            for row in range(start_row + 1, end_row):
                row_start = row * self.row_length
                row_end = (row + 1) * self.row_length - 1
                aa_indices.extend(range(row_start, min(row_end, end) + 1))
            
            # 最后一行
            if end >= end_row * self.row_length:
                aa_indices.extend(range(end_row * self.row_length, end + 1))
        
        # 使用全局segment_index确保跨行连续计数
        for segment_index, aa_index in enumerate(aa_indices):
            self._draw_helix_unit(ax, aa_index, helix_type, segment_index)
    
    def _draw_helix_unit(self, ax, aa_index, helix_type, segment_index):
        """绘制单个螺旋单元"""
        x_center, row = self.aa_positions[aa_index]
        y_center = self.row_y_positions[row]
        
        # 应用位置补偿和间隙
        x_adjustment = -self.position_compensation
        x_adjustment -= segment_index * self.helix_inter_unit_gap
        
        adjusted_x = x_center + x_adjustment
        
        if helix_type == 'H':
            color_left = self.colors['H']  # α-helix颜色
            color_right = self._adjust_color(self.colors['H'], 1.2) # 更亮的颜色
        elif helix_type == 'I':
            color_left = self.colors['I']  # π-helix颜色
            color_right = self._adjust_color(self.colors['I'], 1.2) # 更亮的颜色
        else:  # 3₁₀-helix
            color_left = self.colors['G']  # 3₁₀-helix颜色
            color_right = self._adjust_color(self.colors['G'], 1.2) # 更亮的颜色
        
        # 基本参数
        width = self.parallelogram_width
        overlap = self.helix_intra_unit_overlap
        half_h = self.helix_unit_height / 2
        skew = self.helix_skew
        
        # 左平行四边形（深色）
        left_vertices = [
            (adjusted_x - width, y_center - half_h),
            (adjusted_x - width - skew, y_center + half_h),
            (adjusted_x - skew, y_center + half_h),
            (adjusted_x, y_center - half_h)
        ]
        
        # 右平行四边形（浅色） - 应用重叠
        right_vertices = [
            (adjusted_x - overlap, y_center - half_h),
            (adjusted_x - overlap + skew, y_center + half_h),
            (adjusted_x - overlap + width + skew, y_center + half_h),
            (adjusted_x - overlap + width, y_center - half_h)
        ]
        
        # 绘制顺序：右然后左
        self._draw_parallelogram(ax, right_vertices, color_right, alpha=0.7, zorder=3)
        self._draw_parallelogram(ax, left_vertices, color_left, alpha=0.9, zorder=4)
    
    def _adjust_color(self, color, factor):
        """调整颜色亮度"""
        import colorsys
        # 将十六进制颜色转换为RGB
        r, g, b = [int(color[i:i+2], 16) / 255.0 for i in (1, 3, 5)]
        # 转换为HSV
        h, s, v = colorsys.rgb_to_hsv(r, g, b)
        # 调整亮度
        v = min(1.0, v * factor)
        # 转换回RGB
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        # 转换回十六进制
        return '#{:02x}{:02x}{:02x}'.format(int(r*255), int(g*255), int(b*255))
    
    def _draw_beta_strand(self, ax, start, end):
        """绘制β折叠段"""
        start_x, start_row = self.aa_positions[start]
        end_x, end_row = self.aa_positions[end]
        y_pos = self.row_y_positions[start_row]
        
        length = end - start + 1
        total_length = length * self.aa_spacing
        is_last = self._is_last_beta_segment(end)
        
        if is_last:
            full_length = self._get_beta_total_length(end)
            min_required = self.beta_body_min_length
            
            if full_length * self.aa_spacing < min_required:
                body_length = total_length * 0.7
                head_length = total_length * 0.3
            else:
                head_length = min(
                    self.beta_head_std_length,
                    self.aa_spacing * 1.15,
                    total_length * 0.34,
                )
                body_length = total_length - head_length
        else:
            body_length = total_length
            head_length = 0
        
        body_left_x = start_x - self.aa_spacing / 2
        body_right_x = body_left_x + body_length
        body_upper_y = y_pos + self.beta_body_width / 2
        body_lower_y = y_pos - self.beta_body_width / 2
        
        if is_last and head_length > 0:
            arrow_tip_x = body_right_x + head_length
            arrow_upper_y = y_pos + self.beta_head_width / 2
            arrow_lower_y = y_pos - self.beta_head_width / 2
        else:
            arrow_tip_x = body_right_x
            arrow_upper_y = body_upper_y
            arrow_lower_y = body_lower_y
        
        vertices = [
            (body_left_x, body_lower_y),
            (body_right_x, body_lower_y),
        ]
        
        if is_last and head_length > 0:
            vertices.extend([
                (body_right_x, arrow_lower_y),
                (arrow_tip_x, y_pos),
                (body_right_x, arrow_upper_y),
            ])
        
        vertices.extend([
            (body_right_x, body_upper_y),
            (body_left_x, body_upper_y),
            (body_left_x, body_lower_y)
        ])
        
        codes = [Path.MOVETO] + [Path.LINETO]*(len(vertices)-2) + [Path.CLOSEPOLY]
        arrow_path = Path(vertices, codes)
        
        arrow_patch = patches.PathPatch(
            arrow_path, facecolor=self.colors['E'],
            edgecolor=self.colors['edge'], alpha=0.95, linewidth=0.8, zorder=3
        )
        ax.add_patch(arrow_patch)
    
    def _get_beta_region(self, pos):
        """获取位置的β折叠区域"""
        for (start, end, start_row, end_row) in self.beta_contiguous_regions:
            if start <= pos <= end:
                return (start, end, start_row, end_row)
        return None
    
    def _is_last_beta_segment(self, pos):
        """检查位置是否是β折叠段的最后一个"""
        beta_region = self._get_beta_region(pos)
        return beta_region and pos == beta_region[1]
    
    def _get_beta_total_length(self, pos):
        """获取β折叠段的总长度"""
        beta_region = self._get_beta_region(pos)
        return beta_region[1] - beta_region[0] + 1 if beta_region else 1
    
    def _draw_bend(self, ax, start, end):
        """绘制弯曲"""
        start_x, start_row = self.aa_positions[start]
        end_x, end_row = self.aa_positions[end]
        y_pos = self.row_y_positions[start_row]
        
        line_start_x = start_x - self.aa_spacing / 2
        line_end_x = end_x + self.aa_spacing / 2
        
        ax.plot([line_start_x, line_end_x], [y_pos, y_pos], 
                color=self.colors['S'], linewidth=2.5, solid_capstyle='round', zorder=3)
    
    def _draw_turn(self, ax, start, end):
        """绘制转角"""
        start_x, start_row = self.aa_positions[start]
        end_x, end_row = self.aa_positions[end]
        y_pos = self.row_y_positions[start_row]
        
        line_start_x = start_x - self.aa_spacing / 2
        line_end_x = end_x + self.aa_spacing / 2
        
        x = np.linspace(line_start_x, line_end_x, self.smooth_points)
        amplitude = min(12, 6 + (end - start + 1))
        y = y_pos + amplitude * np.sin(np.pi * (x - line_start_x) / (line_end_x - line_start_x))
        
        ax.plot(x, y, color=self.colors['T'], linewidth=2.5, solid_capstyle='round', zorder=3)
    
    def _draw_coil(self, ax, start, end):
        start_x, start_row = self.aa_positions[start]
        end_x, end_row = self.aa_positions[end]
        y_pos = self.row_y_positions[start_row]
        
        line_start_x = start_x - self.aa_spacing / 2
        line_end_x = end_x + self.aa_spacing / 2
        
        ax.plot([line_start_x, line_end_x], [y_pos, y_pos], 
                color=self.colors['C'], linewidth=2.5, solid_capstyle='round', zorder=3)
    
    def _adjust_layout(self):
        """调整图形布局"""
        plt.subplots_adjust(
            left=self.left_margin,
            right=1.0 - self.right_margin,
            bottom=self.bottom_margin,
            top=1.0 - self.top_margin,
            wspace=self.horizontal_spacing,
            hspace=self.vertical_spacing
        )

    def show(self, fig=None):
        """
        显示图形并确保背景色正确，保证显示与保存一致
        
        Args:
            fig: matplotlib Figure 对象 (可选)，如果未提供则使用当前实例的图形
        """
        # 使用当前实例的图形或传入的图形
        display_fig = fig if fig is not None else self.fig
        
        if display_fig is None:
            raise ValueError("No figure to show. Call create_figure first.")
        
        # 设置图形背景色
        display_fig.set_facecolor(self.colors['background'])
        
        # 设置所有坐标轴背景色
        for ax in display_fig.get_axes():
            ax.set_facecolor(self.colors['background'])
        
        # 应用所有其他视觉设置以确保一致性
        self._apply_consistent_settings(display_fig)
        
        # 显示图形
        plt.show()

    def _apply_consistent_settings(self, fig):
        """应用所有视觉设置以确保显示与保存一致"""
        # 设置所有文本元素的颜色
        for ax in fig.get_axes():
            # 设置标题颜色
            if ax.get_title():
                ax.title.set_color(self.colors['sequence'])
            
            # 设置坐标轴标签颜色
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_color(self.colors['residue_number'])
            
            # 设置图例文本颜色
            legend = ax.get_legend()
            if legend:
                for text in legend.get_texts():
                    text.set_color(self.colors['legend_text'])
            
            # 设置所有文本对象的颜色
            for text in ax.texts:
                # 根据文本内容判断类型
                text_str = text.get_text()
                if text_str.isdigit() and len(text_str) <= 3:  # 残基编号
                    text.set_color(self.colors['residue_number'])
                elif len(text_str) == 1:  # 氨基酸字母
                    text.set_color(self.colors['sequence'])
                elif "Secondary Structure" in text_str:  # 标题
                    text.set_color(self.colors['sequence'])
                else:  # 其他文本
                    text.set_color(self.colors['legend_text'])
        
        # 设置所有补丁对象的边缘颜色
        for ax in fig.get_axes():
            for patch in ax.patches:
                if isinstance(patch, patches.Rectangle) or isinstance(patch, patches.PathPatch):
                    patch.set_edgecolor(self.colors['edge'])
        
        # 设置所有线条对象的颜色
        for ax in fig.get_axes():
            for line in ax.lines:
                # 根据线条位置判断类型
                ydata = line.get_ydata()
                if len(ydata) > 0:
                    y_pos = np.mean(ydata)
                    if abs(y_pos - self.row_y_positions.get(0, 0)) < 5:  # 近似在二级结构位置
                        # 根据颜色判断具体类型
                        color = line.get_color()
                        if color == self.colors['T']:
                            line.set_linewidth(2.5)
                        elif color == self.colors['S']:
                            line.set_linewidth(2.5)
                        elif color == self.colors['C']:
                            line.set_linewidth(2.5)
