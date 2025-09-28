"""
Amino acid properties and characteristics
"""

AA_PROPERTIES = {
    'A': {'name': 'Alanine', 'type': 'hydrophobic', 'weight': 89.1},
    'C': {'name': 'Cysteine', 'type': 'polar', 'weight': 121.2},
    'D': {'name': 'Aspartic Acid', 'type': 'acidic', 'weight': 133.1},
    'E': {'name': 'Glutamic Acid', 'type': 'acidic', 'weight': 147.1},
    'F': {'name': 'Phenylalanine', 'type': 'hydrophobic', 'weight': 165.2},
    'G': {'name': 'Glycine', 'type': 'special', 'weight': 75.1},
    'H': {'name': 'Histidine', 'type': 'basic', 'weight': 155.2},
    'I': {'name': 'Isoleucine', 'type': 'hydrophobic', 'weight': 131.2},
    'K': {'name': 'Lysine', 'type': 'basic', 'weight': 146.2},
    'L': {'name': 'Leucine', 'type': 'hydrophobic', 'weight': 131.2},
    'M': {'name': 'Methionine', 'type': 'hydrophobic', 'weight': 149.2},
    'N': {'name': 'Asparagine', 'type': 'polar', 'weight': 132.1},
    'P': {'name': 'Proline', 'type': 'special', 'weight': 115.1},
    'Q': {'name': 'Glutamine', 'type': 'polar', 'weight': 146.2},
    'R': {'name': 'Arginine', 'type': 'basic', 'weight': 174.2},
    'S': {'name': 'Serine', 'type': 'polar', 'weight': 105.1},
    'T': {'name': 'Threonine', 'type': 'polar', 'weight': 119.1},
    'V': {'name': 'Valine', 'type': 'hydrophobic', 'weight': 117.1},
    'W': {'name': 'Tryptophan', 'type': 'hydrophobic', 'weight': 204.2},
    'Y': {'name': 'Tyrosine', 'type': 'polar', 'weight': 181.2}
}

# Color schemes for amino acids
AA_COLORS = {
    'hydrophobic': '#FF6B6B',
    'polar': '#4ECDC4',
    'acidic': '#FFE66D',
    'basic': '#45B7D1',
    'special': '#96CEB4'
}

# Secondary structure color scheme
SS_COLORS = {
    'H': '#FF6B6B',  # α-helix - red
    'E': '#FFE66D',  # β-strand - yellow
    'C': '#45B7D1',  # coil - blue
    'G': '#FF8E72',  # 3₁₀-helix - orange
    'I': '#C44D58',  # π-helix - dark red
    'T': '#4ECDC4',  # turn - teal
    'S': '#96CEB4'   # bend - green
}