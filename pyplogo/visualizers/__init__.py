def __init__(self, 
             row_length=50,
             amino_acid_spacing=32,
             helix_intra_unit_overlap=12,
             helix_inter_unit_gap=4,
             parallelogram_width=20,
             helix_unit_height=8,
             helix_skew=16,
             beta_body_width_ratio=0.3,
             beta_head_width_ratio=0.45,
             beta_head_length_ratio=0.6,
             aa_font_size=10,
             aa_font_family='Times New Roman'):
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
    """
    self.row_length = row_length
    self.aa_spacing = amino_acid_spacing
    self.aa_font_size = aa_font_size
    self.aa_font_family = aa_font_family
    
    # Helix parameters
    self.helix_intra_unit_overlap = helix_intra_unit_overlap
    self.helix_inter_unit_gap = helix_inter_unit_gap
    self.parallelogram_width = parallelogram_width
    self.helix_unit_height = helix_unit_height
    self.helix_skew = helix_skew
    
    # Beta strand parameters
    self.beta_body_width_ratio = beta_body_width_ratio
    self.beta_head_width_ratio = beta_head_width_ratio
    self.beta_head_length_ratio = beta_head_length_ratio
    
    # 其余初始化代码保持不变...
    self._calculate_derived_parameters()
    self._setup_style_parameters()