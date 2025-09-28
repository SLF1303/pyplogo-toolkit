"""
PyPLogo - Protein Secondary Structure Visualization Tool
A Python package for extracting and visualizing protein secondary structures
from PDB/CIF files using DSSP with publication-quality output.
"""

__version__ = "0.1.0"
__author__ = "Your Name"
__email__ = "your.email@example.com"

# Import main classes for easier access
from .visualizers.secondary_structure import SecondaryStructureVisualizer
from .visualizers.themes import ScientificTheme, NatureTheme, \
    DesertTheme, ArcticTheme, TropicalTheme, GemstoneTheme, VintageTheme, \
    SpaceTheme, SpringTheme, AutumnTheme, CoralReefTheme, AuroraTheme, \
    MetalTheme, CandyTheme, OasisTheme, StarlightTheme, SakuraTheme, \
    VolcanoTheme, JadeTheme, TwilightTheme, OnyxTheme, LavenderTheme, \
    HoneyTheme, MintTheme, CoralTheme, AmethystTheme, LemonTheme, \
    SapphireTheme, RubyTheme, EmeraldTheme, TopazTheme, MalachiteTheme, \
    OpalTheme, PearlTheme

# Define public API
__all__ = [
    'SecondaryStructureVisualizer',
    'ScientificTheme',
    'NatureTheme',
    'DesertTheme',
    'ArcticTheme',
    'TropicalTheme',
    'GemstoneTheme',
    'VintageTheme',
    'SpaceTheme',
    'SpringTheme',
    'AutumnTheme',
    'CoralReefTheme',
    'AuroraTheme',
    'MetalTheme',
    'CandyTheme',
    'OasisTheme',
    'StarlightTheme',
    'SakuraTheme',
    'VolcanoTheme',
    'JadeTheme',
    'TwilightTheme',
    'OnyxTheme',
    'LavenderTheme',
    'HoneyTheme',
    'MintTheme',
    'CoralTheme',
    'AmethystTheme',
    'LemonTheme',
    'SapphireTheme',
    'RubyTheme',
    'EmeraldTheme',
    'TopazTheme',
    'MalachiteTheme',
    'OpalTheme',
    'PearlTheme'
]