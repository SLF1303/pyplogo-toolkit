from dataclasses import dataclass
from typing import Dict, Optional
import matplotlib.colors as mcolors

@dataclass
class Theme:
    """Base theme class for visualization"""
    name: str
    background_color: str
    font_family: str
    font_size: int
    aa_colors: Dict[str, str]
    ss_colors: Dict[str, str]
    ss_styles: Dict[str, Dict]
    sequence_color: str = "#2D3748"
    residue_number_color: str = "#718096"
    disulfide_bond_color: str = "#D32F2F"  # 二硫键颜色
    
    def __post_init__(self):
        """Validate theme parameters"""
        self._validate_colors()
    
    def _validate_colors(self):
        """Validate color codes"""
        for color_name, color_value in self.aa_colors.items():
            if not self._is_valid_color(color_value):
                raise ValueError(f"Invalid color for {color_name}: {color_value}")
        
        for ss_type, color_value in self.ss_colors.items():
            if not self._is_valid_color(color_value):
                raise ValueError(f"Invalid color for {ss_type}: {color_value}")
        
        # 验证新增的颜色属性
        if not self._is_valid_color(self.sequence_color):
            raise ValueError(f"Invalid sequence color: {self.sequence_color}")
        if not self._is_valid_color(self.residue_number_color):
            raise ValueError(f"Invalid residue number color: {self.residue_number_color}")
        if not self._is_valid_color(self.disulfide_bond_color):
            raise ValueError(f"Invalid disulfide bond color: {self.disulfide_bond_color}")
    
    def _is_valid_color(self, color: str) -> bool:
        """Check if color is valid"""
        try:
            mcolors.to_rgba(color)
            return True
        except ValueError:
            return False
    
    def apply_to_visualizer(self, visualizer):
        """将主题设置应用到可视化器"""
        # 设置字体
        visualizer.aa_font_family = self.font_family
        visualizer.aa_font_size = self.font_size
        
        # 创建颜色字典
        colors = {
            'H': self.ss_colors['H'],
            'E': self.ss_colors['E'],
            'G': self.ss_colors['G'],
            'I': self.ss_colors['I'],
            'T': self.ss_colors['T'],
            'S': self.ss_colors['S'],
            'C': self.ss_colors['C'],
            'background': self.background_color,
            'sequence': self.sequence_color,
            'residue_number': self.residue_number_color,
            'disulfide_bond': self.disulfide_bond_color
        }
        
        # 应用颜色
        visualizer.set_colors(colors)

# 科学主题
class ScientificTheme(Theme):
    """科学出版风格主题 - 经典蓝色调"""
    def __init__(self):
        super().__init__(
            name="scientific",
            background_color="#FFFFFF",
            font_family="DejaVu Sans",
            font_size=10,
            aa_colors={
                'hydrophobic': '#1F77B4',
                'polar': '#FF7F0E',
                'acidic': '#2CA02C',
                'basic': '#D62728',
                'special': '#9467BD'
            },
            ss_colors={
                'H': '#1F77B4',  # α-helix - 蓝色
                'E': '#FF7F0E',  # β-strand - 橙色
                'C': '#7F7F7F',  # coil - 灰色
                'G': '#AEC7E8',  # 3₁₀₀-helix - 浅蓝色
                'I': '#8C564B',  # π-helix - 棕色
                'T': '#2CA02C',  # turn - 绿色
                'S': '#9467BD'   # bend - 紫色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#D32F2F"
        )

# 自然主题
class NatureTheme(Theme):
    """Nature期刊风格主题 - 优雅绿色调"""
    def __init__(self):
        super().__init__(
            name="nature",
            background_color="#FFFFFF",
            font_family="Helvetica",
            font_size=9,
            aa_colors={
                'hydrophobic': '#4E79A7',
                'polar': '#F28E2B',
                'acidic': '#59A14F',
                'basic': '#E15759',
                'special': '#B07AA1'
            },
            ss_colors={
                'H': '#59A14F',  # α-helix - 绿色
                'E': '#E15759',  # β-strand - 红色
                'C': '#BAB0AC',  # coil - 浅灰色
                'G': '#76B7B2',  # 3₁₀₀-helix - 青绿色
                'I': '#D37295',  # π-helix - 粉红色
                'T': '#F28E2B',  # turn - 橙色
                'S': '#4E79A7'   # bend - 蓝色
            },
            ss_styles={
                'H': {'style': 'spiral_3d', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow_ribbon', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'smooth_ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral_narrow', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral_wide', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'smooth_curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'gentle_wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#8B0000"
        )

# 柔和主题
class PastelTheme(Theme):
    """柔和色彩主题 - 淡雅色调"""
    def __init__(self):
        super().__init__(
            name="pastel",
            background_color="#F9F7F7",
            font_family="Avenir",
            font_size=10,
            aa_colors={
                'hydrophobic': '#FFB6C1',
                'polar': '#87CEFA',
                'acidic': '#98FB98',
                'basic': '#DDA0DD',
                'special': '#FFD700'
            },
            ss_colors={
                'H': '#FFB6C1',  # α-helix - 粉红色
                'E': '#87CEFA',  # β-strand - 淡蓝色
                'C': '#D3D3D3',  # coil - 浅灰色
                'G': '#FFD700',  # 3₁₀₀-helix - 金色
                'I': '#DDA0DD',  # π-helix - 淡紫色
                'T': '#98FB98',  # turn - 淡绿色
                'S': '#FFA07A'   # bend - 淡橙色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.1, 'height': 0.7},
                'E': {'style': 'arrow', 'width': 0.9, 'height': 0.5},
                'C': {'style': 'ribbon', 'width': 0.5, 'height': 0.2},
                'G': {'style': 'spiral', 'width': 0.7, 'height': 0.4},
                'I': {'style': 'spiral', 'width': 1.3, 'height': 0.8},
                'T': {'style': 'curve', 'width': 0.6, 'height': 0.3},
                'S': {'style': 'wave', 'width': 0.6, 'height': 0.3}
            },
            disulfide_bond_color="#FF69B4"
        )

# 深色主题
class DarkTheme(Theme):
    """深色主题 - 适合夜间使用"""
    def __init__(self):
        super().__init__(
            name="dark",
            background_color="#2C3E50",
            font_family="Roboto",
            font_size=10,
            sequence_color="#ECF0F1",
            residue_number_color="#BDC3C7",
            aa_colors={
                'hydrophobic': '#E74C3C',
                'polar': '#3498DB',
                'acidic': '#2ECC71',
                'basic': '#9B59B6',
                'special': '#F39C12'
            },
            ss_colors={
                'H': '#E74C3C',  # α-helix - 红色
                'E': '#3498DB',  # β-strand - 蓝色
                'C': '#7F8C8D',  # coil - 灰色
                'G': '#9B59B6',  # 3₁₀₀-helix - 紫色
                'I': '#F39C12',  # π-helix - 橙色
                'T': '#2ECC71',  # turn - 绿色
                'S': '#1ABC9C'   # bend - 青绿色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#FF0000"
        )

# 海洋主题
class OceanTheme(Theme):
    """海洋主题 - 蓝色和绿色调"""
    def __init__(self):
        super().__init__(
            name="ocean",
            background_color="#E8F4F8",
            font_family="Arial",
            font_size=10,
            aa_colors={
                'hydrophobic': '#0077B6',
                'polar': '#00B4D8',
                'acidic': '#90E0EF',
                'basic': '#03045E',
                'special': '#CAF0F8'
            },
            ss_colors={
                'H': '#0077B6',  # α-helix - 深蓝色
                'E': '#00B4D8',  # β-strand - 青色
                'C': '#90E0EF',  # coil - 浅蓝色
                'G': '#03045E',  # 3₁₀₀-helix - 深蓝紫色
                'I': '#CAF0F8',  # π-helix - 淡蓝色
                'T': '#48CAE4',  # turn - 亮蓝色
                'S': '#0096C7'   # bend - 中等蓝色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#000080"
        )

# 大地主题
class EarthTheme(Theme):
    """大地主题 - 棕色和绿色调"""
    def __init__(self):
        super().__init__(
            name="earth",
            background_color="#F5F5DC",
            font_family="Georgia",
            font_size=10,
            aa_colors={
                'hydrophobic': '#8B4513',
                'polar': '#228B22',
                'acidic': '#DAA520',
                'basic': '#A0522D',
                'special': '#556B2F'
            },
            ss_colors={
                'H': '#8B4513',  # α-helix - 棕色
                'E': '#228B22',  # β-strand - 绿色
                'C': '#DAA520',  # coil - 金色
                'G': '#A0522D',  # 3₁₀₀-helix - 红棕色
                'I': '#556B2F',  # π-helix - 暗绿色
                'T': '#6B8E23',  # turn - 橄榄绿色
                'S': '#CD853F'   # bend - 黄褐色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#8B0000"
        )

# 日落主题
class SunsetTheme(Theme):
    """日落主题 - 暖色调"""
    def __init__(self):
        super().__init__(
            name="sunset",
            background_color="#FFF5E1",
            font_family="Comic Sans MS",
            font_size=10,
            aa_colors={
                'hydrophobic': '#FF6B6B',
                'polar': '#FFA726',
                'acidic': '#FFD93D',
                'basic': '#6A0DAD',
                'special': '#FF9E80'
            },
            ss_colors={
                'H': '#FF6B6B',  # α-helix - 红色
                'E': '#FFA726',  # β-strand - 橙色
                'C': '#FFD93D',  # coil - 黄色
                'G': '#6A0DAD',  # 3₁₀₀-helix - 紫色
                'I': '#FF9E80',  # π-helix - 淡橙色
                'T': '#FF7F50',  # turn - 珊瑚色
                'S': '#FF6347'   # bend - 番茄红
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#FF4500"
        )

# 霓虹主题
class NeonTheme(Theme):
    """霓虹主题 - 明亮鲜艳"""
    def __init__(self):
        super().__init__(
            name="neon",
            background_color="#000000",
            font_family="Impact",
            font_size=10,
            sequence_color="#FFFFFF",
            residue_number_color="#CCCCCC",
            aa_colors={
                'hydrophobic': '#FF00FF',
                'polar': '#00FFFF',
                'acidic': '#00FF00',
                'basic': '#FF0000',
                'special': '#FFFF00'
            },
            ss_colors={
                'H': '#FF00FF',  # α-helix - 粉红色
                'E': '#00FFFF',  # β-strand - 青色
                'C': '#00FF00',  # coil - 绿色
                'G': '#FF0000',  # 3₁₀₀-helix - 红色
                'I': '#FFFF00',  # π-helix - 黄色
                'T': '#FF9900',  # turn - 橙色
                'S': '#9900FF'   # bend - 紫色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#FF00FF"
        )

# 单色主题
class MonochromeTheme(Theme):
    """单色主题 - 灰度色调"""
    def __init__(self):
        super().__init__(
            name="monochrome",
            background_color="#FFFFFF",
            font_family="Courier New",
            font_size=10,
            aa_colors={
                'hydrophobic': '#333333',
                'polar': '#666666',
                'acidic': '#999999',
                'basic': '#CCCCCC',
                'special': '#555555'
            },
            ss_colors={
                'H': '#000000',  # α-helix - 黑色
                'E': '#333333',  # β-strand - 深灰色
                'C': '#666666',  # coil - 中灰色
                'G': '#999999',  # 3₁₀₀-helix - 浅灰色
                'I': '#BBBBBB',  # π-helix - 淡灰色
                'T': '#CCCCCC',  # turn - 更淡灰色
                'S': '#DDDDDD'   # bend - 最淡灰色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#000000"
        )

# 森林主题
class ForestTheme(Theme):
    """森林主题 - 绿色调"""
    def __init__(self):
        super().__init__(
            name="forest",
            background_color="#F0F7EE",
            font_family="Verdana",
            font_size=10,
            aa_colors={
                'hydrophobic': '#2E5E4E',
                'polar': '#3D8B78',
                'acidic': '#6ABF93',
                'basic': '#9BD9B9',
                'special': '#C4E8D1'
            },
            ss_colors={
                'H': '#2E5E4E',  # α-helix - 深绿色
                'E': '#3D8B78',  # β-strand - 绿色
                'C': '#6ABF93',  # coil - 浅绿色
                'G': '#9BD9B9',  # 3₁₀₀-helix - 淡绿色
                'I': '#C4E8D1',  # π-helix - 极淡绿色
                'T': '#4CAF50',  # turn - 亮绿色
                'S': '#8BC34A'   # bend - 黄绿色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#006400"
        )

# 沙漠主题
class DesertTheme(Theme):
    """沙漠主题 - 沙黄色调"""
    def __init__(self):
        super().__init__(
            name="desert",
            background_color="#FDF5E6",
            font_family="Georgia",
            font_size=10,
            aa_colors={
                'hydrophobic': '#CD853F',
                'polar': '#D2B48C',
                'acidic': '#F4A460',
                'basic': '#8B4513',
                'special': '#DEB887'
            },
            ss_colors={
                'H': '#CD853F',  # α-helix - 黄褐色
                'E': '#D2B48C',  # β-strand - 浅褐色
                'C': '#F4A460',  # coil - 沙棕色
                'G': '#8B4513',  # 3₁₀₀-helix - 深棕色
                'I': '#DEB887',  # π-helix - 淡棕色
                'T': '#A0522D',  # turn - 红棕色
                'S': '#BC8F8F'   # bend - 玫瑰棕色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#8B0000"
        )

# 北极主题
class ArcticTheme(Theme):
    """北极主题 - 冷蓝色调"""
    def __init__(self):
        super().__init__(
            name="arctic",
            background_color="#F0F8FF",
            font_family="Arial",
            font_size=10,
            aa_colors={
                'hydrophobic': '#4682B4',
                'polar': '#87CEEB',
                'acidic': '#B0E0E6',
                'basic': '#1E90FF',
                'special': '#ADD8E6'
            },
            ss_colors={
                'H': '#4682B4',  # α-helix - 钢蓝色
                'E': '#87CEEB',  # β-strand - 天蓝色
                'C': '#B0E0E6',  # coil - 淡蓝色
                'G': '#1E90FF',  # 3₁₀₀-helix - 道奇蓝
                'I': '#ADD8E6',  # π-helix - 浅蓝色
                'T': '#00BFFF',  # turn - 深天蓝
                'S': '#5F9EA0'   # bend - 卡其色蓝
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#000080"
        )

# 热带主题
class TropicalTheme(Theme):
    """热带主题 - 鲜艳色调"""
    def __init__(self):
        super().__init__(
            name="tropical",
            background_color="#FFFACD",
            font_family="Comic Sans MS",
            font_size=10,
            aa_colors={
                'hydrophobic': '#FF4500',
                'polar': '#32CD32',
                'acidic': '#FFD700',
                'basic': '#FF1493',
                'special': '#00CED1'
            },
            ss_colors={
                'H': '#FF4500',  # α-helix - 橙红色
                'E': '#32CD32',  # β-strand - 酸橙绿
                'C': '#FFD700',  # coil - 金色
                'G': '#FF1493',  # 3₁₀₀-helix - 深粉色
                'I': '#00CED1',  # π-helix - 深青色
                'T': '#FF6347',  # turn - 番茄红
                'S': '#7CFC00'   # bend - 草坪绿
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#FF0000"
        )

# 宝石主题
class GemstoneTheme(Theme):
    """宝石主题 - 宝石色调"""
    def __init__(self):
        super().__init__(
            name="gemstone",
            background_color="#F8F8FF",
            font_family="Times New Roman",
            font_size=10,
            aa_colors={
                'hydrophobic': '#4169E1',
                'polar': '#20B2AA',
                'acidic': '#9370DB',
                'basic': '#FF4500',
                'special': '#FFD700'
            },
            ss_colors={
                'H': '#4169E1',  # α-helix - 皇家蓝
                'E': '#20B2AA',  # β-strand - 浅海绿
                'C': '#9370DB',  # coil - 中紫色
                'G': '#FF4500',  # 3₁₀₀-helix - 橙红色
                'I': '#FFD700',  # π-helix - 金色
                'T': '#DA70D6',  # turn - 兰花紫
                'S': '#3CB371'   # bend - 中海绿
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#800080"
        )

# 复古主题
class VintageTheme(Theme):
    """复古主题 - 怀旧色调"""
    def __init__(self):
        super().__init__(
            name="vintage",
            background_color="#FAF0E6",
            font_family="Garamond",
            font_size=10,
            aa_colors={
                'hydrophobic': '#8B4513',
                'polar': '#556B2F',
                'acidic': '#CD853F',
                'basic': '#800000',
                'special': '#6B8E23'
            },
            ss_colors={
                'H': '#8B4513',  # α-helix - 马鞍棕
                'E': '#556B2F',  # β-strand - 暗橄榄绿
                'C': '#CD853F',  # coil - 秘鲁色
                'G': '#800000',  # 3₁₀₀-helix - 栗色
                'I': '#6B8E23',  # π-helix - 橄榄土褐色
                'T': '#A0522D',  # turn - 赭色
                'S': '#B8860B'   # bend - 暗金色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#8B0000"
        )

# 太空主题
class SpaceTheme(Theme):
    """太空主题 - 深空色调"""
    def __init__(self):
        super().__init__(
            name="space",
            background_color="#000033",
            font_family="Arial",
            font_size=10,
            sequence_color="#FFFFFF",
            residue_number_color="#CCCCCC",
            aa_colors={
                'hydrophobic': '#4B0082',
                'polar': '#000080',
                'acidic': '#191970',
                'basic': '#800080',
                'special': '#483D8B'
            },
            ss_colors={
                'H': '#4B0082',  # α-helix - 靛蓝色
                'E': '#000080',  # β-strand - 海军蓝
                'C': '#191970',  # coil - 午夜蓝
                'G': '#800080',  # 3₁₀₀-helix - 紫色
                'I': '#483D8B',  # π-helix - 暗板岩蓝
                'T': '#6A5ACD',  # turn - 板岩蓝
                'S': '#7B68EE'   # bend - 中板岩蓝
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#FF00FF"
        )

# 春天主题
class SpringTheme(Theme):
    """春天主题 - 清新色调"""
    def __init__(self):
        super().__init__(
            name="spring",
            background_color="#F0FFF0",
            font_family="Verdana",
            font_size=10,
            aa_colors={
                'hydrophobic': '#98FB98',
                'polar': '#87CEEB',
                'acidic': '#FFB6C1',
                'basic': '#FFD700',
                'special': '#DDA0DD'
            },
            ss_colors={
                'H': '#98FB98',  # α-helix - 淡绿色
                'E': '#87CEEB',  # β-strand - 天蓝色
                'C': '#FFB6C1',  # coil - 浅粉色
                'G': '#FFD700',  # 3₁₀₀-helix - 金色
                'I': '#DDA0DD',  # π-helix - 梅红色
                'T': '#FFA07A',  # turn - 浅鲑鱼色
                'S': '#BA55D3'   # bend - 中兰花紫
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.1, 'height': 0.7},
                'E': {'style': 'arrow', 'width': 0.9, 'height': 0.5},
                'C': {'style': 'ribbon', 'width': 0.5, 'height': 0.2},
                'G': {'style': 'spiral', 'width': 0.7, 'height': 0.4},
                'I': {'style': 'spiral', 'width': 1.3, 'height': 0.8},
                'T': {'style': 'curve', 'width': 0.6, 'height': 0.3},
                'S': {'style': 'wave', 'width': 0.6, 'height': 0.3}
            },
            disulfide_bond_color="#32CD32"
        )

# 秋天主题
class AutumnTheme(Theme):
    """秋天主题 - 温暖色调"""
    def __init__(self):
        super().__init__(
            name="autumn",
            background_color="#FFF8DC",
            font_family="Georgia",
            font_size=10,
            aa_colors={
                'hydrophobic': '#8B4513',
                'polar': '#D2691E',
                'acidic': '#CD853F',
                'basic': '#A52A2A',
                'special': '#B8860B'
            },
            ss_colors={
                'H': '#8B4513',  # α-helix - 马鞍棕
                'E': '#D2691E',  # β-strand - 巧克力色
                'C': '#CD853F',  # coil - 秘鲁色
                'G': '#A52A2A',  # 3₁₀₀-helix - 棕色
                'I': '#B8860B',  # π-helix - 暗金色
                'T': '#DAA520',  # turn - 金色
                'S': '#BC8F8F'   # bend - 玫瑰棕色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#8B0000"
        )

# 以下继续添加更多主题，总共50个...

# 珊瑚礁主题
class CoralReefTheme(Theme):
    """珊瑚礁主题 - 海洋生物色调"""
    def __init__(self):
        super().__init__(
            name="coral_reef",
            background_color="#E0FFFF",
            font_family="Arial",
            font_size=10,
            aa_colors={
                'hydrophobic': '#FF6B6B',
                'polar': '#48D1CC',
                'acidic': '#FFA07A',
                'basic': '#9370DB',
                'special': '#20B2AA'
            },
            ss_colors={
                'H': '#FF6B6B',  # α-helix - 浅珊瑚色
                'E': '#48D1CC',  # β-strand - 绿松石色
                'C': '#FFA07A',  # coil - 浅鲑鱼色
                'G': '#9370DB',  # 3₁₀₀-helix - 中紫色
                'I': '#20B2AA',  # π-helix - 浅海绿
                'T': '#00CED1',  # turn - 深青色
                'S': '#87CEFA'   # bend - 浅天蓝色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.1, 'height': 0.7},
                'E': {'style': 'arrow', 'width': 0.9, 'height': 0.5},
                'C': {'style': 'ribbon', 'width': 0.5, 'height': 0.2},
                'G': {'style': 'spiral', 'width': 0.7, 'height': 0.4},
                'I': {'style': 'spiral', 'width': 1.3, 'height': 0.8},
                'T': {'style': 'curve', 'width': 0.6, 'height': 0.3},
                'S': {'style': 'wave', 'width': 0.6, 'height': 0.3}
            },
            disulfide_bond_color="#FF4500"
        )

# 极光主题
class AuroraTheme(Theme):
    """极光主题 - 北极光色调"""
    def __init__(self):
        super().__init__(
            name="aurora",
            background_color="#000000",
            font_family="Roboto",
            font_size=10,
            sequence_color="#FFFFFF",
            residue_number_color="#CCCCCC",
            aa_colors={
                'hydrophobic': '#00FF7F',
                'polar': '#00FFFF',
                'acidic': '#9370DB',
                'basic': '#FF00FF',
                'special': '#7CFC00'
            },
            ss_colors={
                'H': '#00FF7F',  # α-helix - 春绿色
                'E': '#00FFFF',  # β-strand - 青色
                'C': '#9370DB',  # coil - 中紫色
                'G': '#FF00FF',  # 3₁₀₀-helix - 洋红色
                'I': '#7CFC00',  # π-helix - 草坪绿
                'T': '#40E0D0',  # turn - 绿松石色
                'S': '#BA55D3'   # bend - 中兰花紫
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#FF00FF"
        )

# 金属主题
class MetalTheme(Theme):
    """金属主题 - 工业金属色调"""
    def __init__(self):
        super().__init__(
            name="metal",
            background_color="#D3D3D3",
            font_family="Arial Black",
            font_size=10,
            aa_colors={
                'hydrophobic': '#808080',
                'polar': '#A9A9A9',
                'acidic': '#C0C0C0',
                'basic': '#696969',
                'special': '#2F4F4F'
            },
            ss_colors={
                'H': '#808080',  # α-helix - 灰色
                'E': '#A9A9A9',  # β-strand - 暗灰色
                'C': '#C0C0C0',  # coil - 银色
                'G': '#696969',  # 3₁₀₀-helix - 暗灰色
                'I': '#2F4F4F',  # π-helix - 深石板灰
                'T': '#708090',  # turn - 石板灰
                'S': '#778899'   # bend - 浅石板灰
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#000000"
        )

# 糖果主题
class CandyTheme(Theme):
    """糖果主题 - 甜美色调"""
    def __init__(self):
        super().__init__(
            name="candy",
            background_color="#FFF5EE",
            font_family="Comic Sans MS",
            font_size=10,
            aa_colors={
                'hydrophobic': '#FFB6C1',
                'polar': '#87CEFA',
                'acidic': '#98FB98',
                'basic': '#DDA0DD',
                'special': '#FFD700'
            },
            ss_colors={
                'H': '#FFB6C1',  # α-helix - 浅粉色
                'E': '#87CEFA',  # β-strand - 浅蓝色
                'C': '#98FB98',  # coil - 浅绿色
                'G': '#DDA0DD',  # 3₁₀₀-helix - 浅紫色
                'I': '#FFD700',  # π-helix - 金色
                'T': '#FFA07A',  # turn - 浅鲑鱼色
                'S': '#FF69B4'   # bend - 热粉色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.1, 'height': 0.7},
                'E': {'style': 'arrow', 'width': 0.9, 'height': 0.5},
                'C': {'style': 'ribbon', 'width': 0.5, 'height': 0.2},
                'G': {'style': 'spiral', 'width': 0.7, 'height': 0.4},
                'I': {'style': 'spiral', 'width': 1.3, 'height': 0.8},
                'T': {'style': 'curve', 'width': 0.6, 'height': 0.3},
                'S': {'style': 'wave', 'width': 0.6, 'height': 0.3}
            },
            disulfide_bond_color="#FF1493"
        )

# 沙漠绿洲主题
class OasisTheme(Theme):
    """沙漠绿洲主题 - 蓝绿色调"""
    def __init__(self):
        super().__init__(
            name="oasis",
            background_color="#F0FFF0",
            font_family="Arial",
            font_size=10,
            aa_colors={
                'hydrophobic': '#00CED1',
                'polar': '#20B2AA',
                'acidic': '#AFEEEE',
                'basic': '#008080',
                'special': '#48D1CC'
            },
            ss_colors={
                'H': '#00CED1',  # α-helix - 深青色
                'E': '#20B2AA',  # β-strand - 浅海绿
                'C': '#AFEEEE',  # coil - 苍蓝色
                'G': '#008080',  # 3₁₀₀-helix - 蓝绿色
                'I': '#48D1CC',  # π-helix - 中绿松石色
                'T': '#40E0D0',  # turn - 绿松石色
                'S': '#00FA9A'   # bend - 中春绿色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.1, 'height': 0.7},
                'E': {'style': 'arrow', 'width': 0.9, 'height': 0.5},
                'C': {'style': 'ribbon', 'width': 0.5, 'height': 0.2},
                'G': {'style': 'spiral', 'width': 0.7, 'height': 0.4},
                'I': {'style': 'spiral', 'width': 1.3, 'height': 0.8},
                'T': {'style': 'curve', 'width': 0.6, 'height': 0.3},
                'S': {'style': 'wave', 'width': 0.6, 'height': 0.3}
            },
            disulfide_bond_color="#008B8B"
        )

# 星空主题
class StarlightTheme(Theme):
    """星空主题 - 深蓝与金色调"""
    def __init__(self):
        super().__init__(
            name="starlight",
            background_color="#000033",
            font_family="Arial",
            font_size=10,
            sequence_color="#FFFFFF",
            residue_number_color="#CCCCCC",
            aa_colors={
                'hydrophobic': '#4169E1',
                'polar': '#9370DB',
                'acidic': '#F0E68C',
                'basic': '#FFD700',
                'special': '#87CEEB'
            },
            ss_colors={
                'H': '#4169E1',  # α-helix - 皇家蓝
                'E': '#9370DB',  # β-strand - 中紫色
                'C': '#F0E68C',  # coil - 卡其色
                'G': '#FFD700',  # 3₁₀₀-helix - 金色
                'I': '#87CEEB',  # π-helix - 天蓝色
                'T': '#FFA500',  # turn - 橙色
                'S': '#ADD8E6'   # bend - 浅蓝色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#FFFF00"
        )

# 樱花主题
class SakuraTheme(Theme):
    """樱花主题 - 粉红色调"""
    def __init__(self):
        super().__init__(
            name="sakura",
            background_color="#FFF0F5",
            font_family="MS Gothic",
            font_size=10,
            aa_colors={
                'hydrophobic': '#FF69B4',
                'polar': '#FFB6C1',
                'acidic': '#FFC0CB',
                'basic': '#DB7093',
                'special': '#D8BFD8'
            },
            ss_colors={
                'H': '#FF69B4',  # α-helix - 热粉色
                'E': '#FFB6C1',  # β-strand - 浅粉色
                'C': '#FFC0CB',  # coil - 粉红色
                'G': '#DB7093',  # 3₁₀₀-helix - 苍紫红色
                'I': '#D8BFD8',  # π-helix - 蓟色
                'T': '#FF1493',  # turn - 深粉色
                'S': '#DA70D6'   # bend - 兰花紫
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.1, 'height': 0.7},
                'E': {'style': 'arrow', 'width': 0.9, 'height': 0.5},
                'C': {'style': 'ribbon', 'width': 0.5, 'height': 0.2},
                'G': {'style': 'spiral', 'width': 0.7, 'height': 0.4},
                'I': {'style': 'spiral', 'width': 1.3, 'height': 0.8},
                'T': {'style': 'curve', 'width': 0.6, 'height': 0.3},
                'S': {'style': 'wave', 'width': 0.6, 'height': 0.3}
            },
            disulfide_bond_color="#C71585"
        )

# 火山主题
class VolcanoTheme(Theme):
    """火山主题 - 红黑色调"""
    def __init__(self):
        super().__init__(
            name="volcano",
            background_color="#2F4F4F",
            font_family="Arial Black",
            font_size=10,
            sequence_color="#FFFFFF",
            residue_number_color="#CCCCCC",
            aa_colors={
                'hydrophobic': '#8B0000',
                'polar': '#B22222',
                'acidic': '#DC143C',
                'basic': '#FF4500',
                'special': '#FF6347'
            },
            ss_colors={
                'H': '#8B0000',  # α-helix - 深红色
                'E': '#B22222',  # β-strand - 火砖色
                'C': '#DC143C',  # coil - 深红色
                'G': '#FF4500',  # 3₁₀₀-helix - 橙红色
                'I': '#FF6347',  # π-helix - 番茄红
                'T': '#CD5C5C',  # turn - 印度红
                'S': '#FF0000'   # bend - 红色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#FF4500"
        )

# 翡翠主题
class JadeTheme(Theme):
    """翡翠主题 - 绿色调"""
    def __init__(self):
        super().__init__(
            name="jade",
            background_color="#F0FFF0",
            font_family="Georgia",
            font_size=10,
            aa_colors={
                'hydrophobic': '#2E8B57',
                'polar': '#3CB371',
                'acidic': '#66CDAA',
                'basic': '#8FBC8F',
                'special': '#98FB98'
            },
            ss_colors={
                'H': '#2E8B57',  # α-helix - 海绿色
                'E': '#3CB371',  # β-strand - 中海绿色
                'C': '#66CDAA',  # coil - 中蓝绿色
                'G': '#8FBC8F',  # 3₁₀₀-helix - 暗海绿色
                'I': '#98FB98',  # π-helix - 淡绿色
                'T': '#00FA9A',  # turn - 中春绿色
                'S': '#7CFC00'   # bend - 草坪绿
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.1, 'height': 0.7},
                'E': {'style': 'arrow', 'width': 0.9, 'height': 0.5},
                'C': {'style': 'ribbon', 'width': 0.5, 'height': 0.2},
                'G': {'style': 'spiral', 'width': 0.7, 'height': 0.4},
                'I': {'style': 'spiral', 'width': 1.3, 'height': 0.8},
                'T': {'style': 'curve', 'width': 0.6, 'height': 0.3},
                'S': {'style': 'wave', 'width': 0.6, 'height': 0.3}
            },
            disulfide_bond_color="#006400"
        )

# 黄昏主题
class TwilightTheme(Theme):
    """黄昏主题 - 紫蓝色调"""
    def __init__(self):
        super().__init__(
            name="twilight",
            background_color="#E6E6FA",
            font_family="Times New Roman",
            font_size=10,
            aa_colors={
                'hydrophobic': '#483D8B',
                'polar': '#6A5ACD',
                'acidic': '#9370DB',
                'basic': '#8A2BE2',
                'special': '#7B68EE'
            },
            ss_colors={
                'H': '#483D8B',  # α-helix - 暗板岩蓝
                'E': '#6A5ACD',  # β-strand - 板岩蓝
                'C': '#9370DB',  # coil - 中紫色
                'G': '#8A2BE2',  # 3₁₀₀-helix - 蓝紫色
                'I': '#7B68EE',  # π-helix - 中板岩蓝
                'T': '#9370DB',  # turn - 中紫色
                'S': '#BA55D3'   # bend - 中兰花紫
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#4B0082"
        )

# 玛瑙主题
class OnyxTheme(Theme):
    """玛瑙主题 - 黑灰色调"""
    def __init__(self):
        super().__init__(
            name="onyx",
            background_color="#1C1C1C",
            font_family="Verdana",
            font_size=10,
            sequence_color="#FFFFFF",
            residue_number_color="#CCCCCC",
            aa_colors={
                'hydrophobic': '#696969',
                'polar': '#808080',
                'acidic': '#A9A9A9',
                'basic': '#C0C0C0',
                'special': '#D3D3D3'
            },
            ss_colors={
                'H': '#696969',  # α-helix - 暗灰色
                'E': '#808080',  # β-strand - 灰色
                'C': '#A9A9A9',  # coil - 暗灰色
                'G': '#C0C0C0',  # 3₁₀₀-helix - 银色
                'I': '#D3D3D3',  # π-helix - 浅灰色
                'T': '#BEBEBE',  # turn - 中灰色
                'S': '#DCDCDC'   # bend - 淡灰色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#FFFFFF"
        )

# 薰衣草主题
class LavenderTheme(Theme):
    """薰衣草主题 - 紫色调"""
    def __init__(self):
        super().__init__(
            name="lavender",
            background_color="#F5F0FF",
            font_family="Arial",
            font_size=10,
            aa_colors={
                'hydrophobic': '#9370DB',
                'polar': '#BA55D3',
                'acidic': '#D8BFD8',
                'basic': '#DA70D6',
                'special': '#EE82EE'
            },
            ss_colors={
                'H': '#9370DB',  # α-helix - 中紫色
                'E': '#BA55D3',  # β-strand - 中兰花紫
                'C': '#D8BFD8',  # coil - 蓟色
                'G': '#DA70D6',  # 3₁₀₀-helix - 兰花紫
                'I': '#EE82EE',  # π-helix - 紫罗兰色
                'T': '#DDA0DD',  # turn - 梅红色
                'S': '#E6E6FA'   # bend - 淡紫色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.1, 'height': 0.7},
                'E': {'style': 'arrow', 'width': 0.9, 'height': 0.5},
                'C': {'style': 'ribbon', 'width': 0.5, 'height': 0.2},
                'G': {'style': 'spiral', 'width': 0.7, 'height': 0.4},
                'I': {'style': 'spiral', 'width': 1.3, 'height': 0.8},
                'T': {'style': 'curve', 'width': 0.6, 'height': 0.3},
                'S': {'style': 'wave', 'width': 0.6, 'height': 0.3}
            },
            disulfide_bond_color="#800080"
        )

# 蜂蜜主题
class HoneyTheme(Theme):
    """蜂蜜主题 - 金黄色调"""
    def __init__(self):
        super().__init__(
            name="honey",
            background_color="#FFF8DC",
            font_family="Georgia",
            font_size=10,
            aa_colors={
                'hydrophobic': '#DAA520',
                'polar': '#FFD700',
                'acidic': '#F0E68C',
                'basic': '#B8860B',
                'special': '#FFA500'
            },
            ss_colors={
                'H': '#DAA520',  # α-helix - 金色
                'E': '#FFD700',  # β-strand - 金色
                'C': '#F0E68C',  # coil - 卡其色
                'G': '#B8860B',  # 3₁₀₀-helix - 暗金色
                'I': '#FFA500',  # π-helix - 橙色
                'T': '#CD853F',  # turn - 秘鲁色
                'S': '#D2691E'   # bend - 巧克力色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#8B4513"
        )

# 薄荷主题
class MintTheme(Theme):
    """薄荷主题 - 清新绿色调"""
    def __init__(self):
        super().__init__(
            name="mint",
            background_color="#F5FFFA",
            font_family="Arial",
            font_size=10,
            aa_colors={
                'hydrophobic': '#3CB371',
                'polar': '#2E8B57',
                'acidic': '#98FB98',
                'basic': '#00FA9A',
                'special': '#7FFF00'
            },
            ss_colors={
                'H': '#3CB371',  # α-helix - 中海绿色
                'E': '#2E8B57',  # β-strand - 海绿色
                'C': '#98FB98',  # coil - 淡绿色
                'G': '#00FA9A',  # 3₁₀₀-helix - 中春绿色
                'I': '#7FFF00',  # π-helix - 查特酒绿
                'T': '#00FF7F',  # turn - 春绿色
                'S': '#7CFC00'   # bend - 草坪绿
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.1, 'height': 0.7},
                'E': {'style': 'arrow', 'width': 0.9, 'height': 0.5},
                'C': {'style': 'ribbon', 'width': 0.5, 'height': 0.2},
                'G': {'style': 'spiral', 'width': 0.7, 'height': 0.4},
                'I': {'style': 'spiral', 'width': 1.3, 'height': 0.8},
                'T': {'style': 'curve', 'width': 0.6, 'height': 0.3},
                'S': {'style': 'wave', 'width': 0.6, 'height': 0.3}
            },
            disulfide_bond_color="#008000"
        )

# 珊瑚主题
class CoralTheme(Theme):
    """珊瑚主题 - 珊瑚色调"""
    def __init__(self):
        super().__init__(
            name="coral",
            background_color="#FFF5EE",
            font_family="Arial",
            font_size=10,
            aa_colors={
                'hydrophobic': '#FF6347',
                'polar': '#FF7F50',
                'acidic': '#FA8072',
                'basic': '#F08080',
                'special': '#CD5C5C'
            },
            ss_colors={
                'H': '#FF6347',  # α-helix - 番茄红
                'E': '#FF7F50',  # β-strand - 珊瑚色
                'C': '#FA8072',  # coil - 鲑鱼色
                'G': '#F08080',  # 3₁₀₀-helix - 浅珊瑚色
                'I': '#CD5C5C',  # π-helix - 印度红
                'T': '#E9967A',  # turn - 深鲑鱼色
                'S': '#FFA07A'   # bend - 浅鲑鱼色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.1, 'height': 0.7},
                'E': {'style': 'arrow', 'width': 0.9, 'height': 0.5},
                'C': {'style': 'ribbon', 'width': 0.5, 'height': 0.2},
                'G': {'style': 'spiral', 'width': 0.7, 'height': 0.4},
                'I': {'style': 'spiral', 'width': 1.3, 'height': 0.8},
                'T': {'style': 'curve', 'width': 0.6, 'height': 0.3},
                'S': {'style': 'wave', 'width': 0.6, 'height': 0.3}
            },
            disulfide_bond_color="#B22222"
        )

# 紫水晶主题
class AmethystTheme(Theme):
    """紫水晶主题 - 紫色调"""
    def __init__(self):
        super().__init__(
            name="amethyst",
            background_color="#F5F0FF",
            font_family="Times New Roman",
            font_size=10,
            aa_colors={
                'hydrophobic': '#9370DB',
                'polar': '#8A2BE2',
                'acidic': '#9932CC',
                'basic': '#BA55D3',
                'special': '#DA70D6'
            },
            ss_colors={
                'H': '#9370DB',  # α-helix - 中紫色
                'E': '#8A2BE2',  # β-strand - 蓝紫色
                'C': '#9932CC',  # coil - 深兰花紫
                'G': '#BA55D3',  # 3₁₀₀-helix - 中兰花紫
                'I': '#DA70D6',  # π-helix - 兰花紫
                'T': '#D8BFD8',  # turn - 蓟色
                'S': '#EE82EE'   # bend - 紫罗兰色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#4B0082"
        )

# 柠檬主题
class LemonTheme(Theme):
    """柠檬主题 - 黄绿色调"""
    def __init__(self):
        super().__init__(
            name="lemon",
            background_color="#FFFFF0",
            font_family="Arial",
            font_size=10,
            aa_colors={
                'hydrophobic': '#9ACD32',
                'polar': '#ADFF2F',
                'acidic': '#7FFF00',
                'basic': '#32CD32',
                'special': '#00FF00'
            },
            ss_colors={
                'H': '#9ACD32',  # α-helix - 黄绿色
                'E': '#ADFF2F',  # β-strand - 绿黄色
                'C': '#7FFF00',  # coil - 查特酒绿
                'G': '#32CD32',  # 3₁₀₀-helix - 酸橙绿
                'I': '#00FF00',  # π-helix - 绿色
                'T': '#7CFC00',  # turn - 草坪绿
                'S': '#00FF7F'   # bend - 春绿色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.1, 'height': 0.7},
                'E': {'style': 'arrow', 'width': 0.9, 'height': 0.5},
                'C': {'style': 'ribbon', 'width': 0.5, 'height': 0.2},
                'G': {'style': 'spiral', 'width': 0.7, 'height': 0.4},
                'I': {'style': 'spiral', 'width': 1.3, 'height': 0.8},
                'T': {'style': 'curve', 'width': 0.6, 'height': 0.3},
                'S': {'style': 'wave', 'width': 0.6, 'height': 0.3}
            },
            disulfide_bond_color="#008000"
        )

# 蓝宝石主题
class SapphireTheme(Theme):
    """蓝宝石主题 - 深蓝色调"""
    def __init__(self):
        super().__init__(
            name="sapphire",
            background_color="#F0F8FF",
            font_family="Arial",
            font_size=10,
            aa_colors={
                'hydrophobic': '#000080',
                'polar': '#0000CD',
                'acidic': '#4169E1',
                'basic': '#1E90FF',
                'special': '#4682B4'
            },
            ss_colors={
                'H': '#000080',  # α-helix - 海军蓝
                'E': '#0000CD',  # β-strand - 中蓝色
                'C': '#4169E1',  # coil - 皇家蓝
                'G': '#1E90FF',  # 3₁₀₀-helix - 道奇蓝
                'I': '#4682B4',  # π-helix - 钢蓝色
                'T': '#6495ED',  # turn - 矢车菊蓝
                'S': '#87CEFA'   # bend - 浅天蓝色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#0000FF"
        )

# 红宝石主题
class RubyTheme(Theme):
    """红宝石主题 - 深红色调"""
    def __init__(self):
        super().__init__(
            name="ruby",
            background_color="#FFF0F5",
            font_family="Georgia",
            font_size=10,
            aa_colors={
                'hydrophobic': '#8B0000',
                'polar': '#B22222',
                'acidic': '#DC143C',
                'basic': '#FF0000',
                'special': '#CD5C5C'
            },
            ss_colors={
                'H': '#8B0000',  # α-helix - 深红色
                'E': '#B22222',  # β-strand - 火砖色
                'C': '#DC143C',  # coil - 深红色
                'G': '#FF0000',  # 3₁₀₀-helix - 红色
                'I': '#CD5C5C',  # π-helix - 印度红
                'T': '#FF4500',  # turn - 橙红色
                'S': '#FF6347'   # bend - 番茄红
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#B22222"
        )

# 翡翠主题
class EmeraldTheme(Theme):
    """翡翠主题 - 绿色调"""
    def __init__(self):
        super().__init__(
            name="emerald",
            background_color="#F0FFF0",
            font_family="Verdana",
            font_size=10,
            aa_colors={
                'hydrophobic': '#006400',
                'polar': '#228B22',
                'acidic': '#2E8B57',
                'basic': '#3CB371',
                'special': '#00FA9A'
            },
            ss_colors={
                'H': '#006400',  # α-helix - 深绿色
                'E': '#228B22',  # β-strand - 森林绿
                'C': '#2E8B57',  # coil - 海绿色
                'G': '#3CB371',  # 3₁₀₀-helix - 中海绿色
                'I': '#00FA9A',  # π-helix - 中春绿色
                'T': '#32CD32',  # turn - 酸橙绿
                'S': '#7CFC00'   # bend - 草坪绿
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.1, 'height': 0.7},
                'E': {'style': 'arrow', 'width': 0.9, 'height': 0.5},
                'C': {'style': 'ribbon', 'width': 0.5, 'height': 0.2},
                'G': {'style': 'spiral', 'width': 0.7, 'height': 0.4},
                'I': {'style': 'spiral', 'width': 1.3, 'height': 0.8},
                'T': {'style': 'curve', 'width': 0.6, 'height': 0.3},
                'S': {'style': 'wave', 'width': 0.6, 'height': 0.3}
            },
            disulfide_bond_color="#008000"
        )

# 黄玉主题
class TopazTheme(Theme):
    """黄玉主题 - 黄褐色调"""
    def __init__(self):
        super().__init__(
            name="topaz",
            background_color="#FFF8DC",
            font_family="Georgia",
            font_size=10,
            aa_colors={
                'hydrophobic': '#D2691E',
                'polar': '#CD853F',
                'acidic': '#DAA520',
                'basic': '#B8860B',
                'special': '#F0E68C'
            },
            ss_colors={
                'H': '#D2691E',  # α-helix - 巧克力色
                'E': '#CD853F',  # β-strand - 秘鲁色
                'C': '#DAA520',  # coil - 金色
                'G': '#B8860B',  # 3₁₀₀-helix - 暗金色
                'I': '#F0E68C',  # π-helix - 卡其色
                'T': '#D2B48C',  # turn - 黄褐色
                'S': '#DEB887'   # bend - 浅黄褐色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#8B4513"
        )

# 孔雀石主题
class MalachiteTheme(Theme):
    """孔雀石主题 - 绿色调"""
    def __init__(self):
        super().__init__(
            name="malachite",
            background_color="#F0FFF0",
            font_family="Arial",
            font_size=10,
            aa_colors={
                'hydrophobic': '#008000',
                'polar': '#00FF00',
                'acidic': '#32CD32',
                'basic': '#00FA9A',
                'special': '#7CFC00'
            },
            ss_colors={
                'H': '#008000',  # α-helix - 绿色
                'E': '#00FF00',  # β-strand - 绿色
                'C': '#32CD32',  # coil - 酸橙绿
                'G': '#00FA9A',  # 3₁₀₀-helix - 中春绿色
                'I': '#7CFC00',  # π-helix - 草坪绿
                'T': '#00FF7F',  # turn - 春绿色
                'S': '#ADFF2F'   # bend - 绿黄色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.1, 'height': 0.7},
                'E': {'style': 'arrow', 'width': 0.9, 'height': 0.5},
                'C': {'style': 'ribbon', 'width': 0.5, 'height': 0.2},
                'G': {'style': 'spiral', 'width': 0.7, 'height': 0.4},
                'I': {'style': 'spiral', 'width': 1.3, 'height': 0.8},
                'T': {'style': 'curve', 'width': 0.6, 'height': 0.3},
                'S': {'style': 'wave', 'width': 0.6, 'height': 0.3}
            },
            disulfide_bond_color="#006400"
        )

# 蛋白石主题
class OpalTheme(Theme):
    """蛋白石主题 - 彩虹色调"""
    def __init__(self):
        super().__init__(
            name="opal",
            background_color="#F0F8FF",
            font_family="Arial",
            font_size=10,
            aa_colors={
                'hydrophobic': '#FF4500',
                'polar': '#FFD700',
                'acidic': '#32CD32',
                'basic': '#1E90FF',
                'special': '#9370DB'
            },
            ss_colors={
                'H': '#FF4500',  # α-helix - 橙红色
                'E': '#FFD700',  # β-strand - 金色
                'C': '#32CD32',  # coil - 酸橙绿
                'G': '#1E90FF',  # 3₁₀₀-helix - 道奇蓝
                'I': '#9370DB',  # π-helix - 中紫色
                'T': '#FF69B4',  # turn - 热粉色
                'S': '#00CED1'   # bend - 深青色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.2, 'height': 0.8},
                'E': {'style': 'arrow', 'width': 1.0, 'height': 0.6},
                'C': {'style': 'ribbon', 'width': 0.6, 'height': 0.3},
                'G': {'style': 'spiral', 'width': 0.8, 'height': 0.5},
                'I': {'style': 'spiral', 'width': 1.4, 'height': 0.9},
                'T': {'style': 'curve', 'width': 0.7, 'height': 0.4},
                'S': {'style': 'wave', 'width': 0.7, 'height': 0.4}
            },
            disulfide_bond_color="#FF0000"
        )

# 珍珠主题
class PearlTheme(Theme):
    """珍珠主题 - 柔和色调"""
    def __init__(self):
        super().__init__(
            name="pearl",
            background_color="#FAF0E6",
            font_family="Georgia",
            font_size=10,
            aa_colors={
                'hydrophobic': '#D3D3D3',
                'polar': '#E6E6FA',
                'acidic': '#F5F5DC',
                'basic': '#FFF0F5',
                'special': '#F0F8FF'
            },
            ss_colors={
                'H': '#D3D3D3',  # α-helix - 浅灰色
                'E': '#E6E6FA',  # β-strand - 淡紫色
                'C': '#F5F5DC',  # coil - 浅米色
                'G': '#FFF0F5',  # 3₁₀₀-helix - 淡粉色
                'I': '#F0F8FF',  # π-helix - 淡蓝色
                'T': '#F5F5F5',  # turn - 白色
                'S': '#F0FFF0'   # bend - 淡绿色
            },
            ss_styles={
                'H': {'style': 'spiral', 'width': 1.1, 'height': 0.7},
                'E': {'style': 'arrow', 'width': 0.9, 'height': 0.5},
                'C': {'style': 'ribbon', 'width': 0.5, 'height': 0.2},
                'G': {'style': 'spiral', 'width': 0.7, 'height': 0.4},
                'I': {'style': 'spiral', 'width': 1.3, 'height': 0.8},
                'T': {'style': 'curve', 'width': 0.6, 'height': 0.3},
                'S': {'style': 'wave', 'width': 0.6, 'height': 0.3}
            },
            disulfide_bond_color="#C0C0C0"
        )

# 获取所有主题
def get_all_themes():
    """获取所有主题"""
    return [
        "scientific",
        "nature",
        "pastel",
        "dark",
        "ocean",
        "earth",
        "sunset",
        "neon",
        "monochrome",
        "forest",
        "desert",
        "arctic",
        "tropical", 
        "gemstone",
        "vintage",
        "space",
        "spring",
        "autumn",
        "coral_reef",
        "aurora",
        "metal",
        "candy",
        "oasis",
        "starlight",
        "sakura",
        "volcano",
        "jade",
        "twilight",
        "onyx",
        "lavender",
        "honey",
        "mint",
        "coral",
        "amethyst",
        "lemon",
        "sapphire",
        "ruby",
        "emerald",
        "topaz",
        "malachite",
        "opal",
        "pearl"
    ]

# 主题名称映射字典
THEME_MAP = {
    "scientific": ScientificTheme,
    "nature": NatureTheme,
    "pastel": PastelTheme,
    "dark": DarkTheme,
    "ocean": OceanTheme,
    "earth": EarthTheme,
    "sunset": SunsetTheme,
    "neon": NeonTheme,
    "monochrome": MonochromeTheme,
    "forest": ForestTheme,
    "desert": DesertTheme,
    "arctic": ArcticTheme,
    "tropical": TropicalTheme,
    "gemstone": GemstoneTheme,
    "vintage": VintageTheme,
    "space": SpaceTheme,
    "spring": SpringTheme,
    "autumn": AutumnTheme,
    "coral_reef": CoralReefTheme,
    "aurora": AuroraTheme,
    "metal": MetalTheme,
    "candy": CandyTheme,
    "oasis": OasisTheme,
    "starlight": StarlightTheme,
    "sakura": SakuraTheme,
    "volcano": VolcanoTheme,
    "jade": JadeTheme,
    "twilight": TwilightTheme,
    "onyx": OnyxTheme,
    "lavender": LavenderTheme,
    "honey": HoneyTheme,
    "mint": MintTheme,
    "coral": CoralTheme,
    "amethyst": AmethystTheme,
    "lemon": LemonTheme,
    "sapphire": SapphireTheme,
    "ruby": RubyTheme,
    "emerald": EmeraldTheme,
    "topaz": TopazTheme,
    "malachite": MalachiteTheme,
    "opal": OpalTheme,
    "pearl": PearlTheme
}

def get_theme(name, **overrides):
    """按名称获取主题并应用覆盖设置"""
    if name not in THEME_MAP:
        print(f"警告：未找到主题 '{name}'，使用默认科学主题")
        theme_class = ScientificTheme
    else:
        theme_class = THEME_MAP[name]
    
    # 创建主题实例
    theme = theme_class()
    
    # 应用覆盖设置
    for key, value in overrides.items():
        if hasattr(theme, key):
            setattr(theme, key, value)
        else:
            print(f"警告：主题没有属性 '{key}'，跳过设置")
    
    return theme