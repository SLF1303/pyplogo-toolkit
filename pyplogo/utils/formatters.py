"""
Utility functions for formatting and formatting protein sequences and structures.
"""

from typing import List, Tuple, Dict, Union
import numpy as np

def format_sequence(sequence: str, per_line: int = 50) -> List[str]:
    """
    Format amino acid sequence into multiple lines
    
    Args:
        sequence: Amino acid sequence
        per_line: Number of amino acids per line
        
    Returns:
        List of formatted sequence lines
    """
    if not sequence:
        return []
    
    # 将序列分割成指定长度的块
    chunks = []
    for i in range(0, len(sequence), per_line):
        chunk = sequence[i:i + per_line]
        chunks.append(chunk)
    
    return chunks

def format_secondary_structure(ss_string: str, line_length: int = 50) -> List[str]:
    """
    Format secondary structure string into multiple lines
    
    Args:
        ss_string: Secondary structure string
        line_length: Number of characters per line
        
    Returns:
        List of formatted secondary structure lines
    """
    return [ss_string[i:i+line_length] 
            for i in range(0, len(ss_string), line_length)]

def format_confidence_scores(confidence: np.ndarray, 
                            line_length: int = 50) -> List[np.ndarray]:
    """
    Format confidence scores into multiple lines
    
    Args:
        confidence: Confidence scores array
        line_length: Number of scores per line
        
    Returns:
        List of formatted confidence arrays
    """
    return [confidence[i:i+line_length] 
            for i in range(0, len(confidence), line_length)]

def create_sequence_chunks(sequence: str, chunk_size: int = 50) -> List[Tuple[int, int]]:
    """
    Create chunks for sequence processing
    
    Args:
        sequence: Amino acid sequence
        chunk_size: Size of each chunk
        
    Returns:
        List of (start, end) tuples for each chunk
    """
    return [(i, min(i + chunk_size, len(sequence))) 
            for i in range(0, len(sequence), chunk_size)]

def calculate_sequence_position(residue_index: int, 
                              chunk_size: int = 50, 
                              line_index: int = 0) -> Tuple[int, int]:
    """
    Calculate visual position for a residue
    
    Args:
        residue_index: Index of residue in sequence
        chunk_size: Number of residues per line
        line_index: Starting line index
        
    Returns:
        Tuple of (line_number, position_in_line)
    """
    line_number = line_index + (residue_index // chunk_size)
    position_in_line = residue_index % chunk_size
    return line_number, position_in_line

def format_legend(ss_types: List[str], 
                color_map: Dict[str, str]) -> str:
    """
    Create formatted legend for secondary structure types
    
    Args:
        ss_types: List of secondary structure types
        color_map: Mapping from SS type to color code
        
    Returns:
        Formatted legend string
    """
    legend_lines = []
    for ss_type in ss_types:
        color = color_map.get(ss_type, '#000000')
        legend_lines.append(f"{ss_type}: {color}")
    return "\n".join(legend_lines)

def normalize_confidence_scores(scores: np.ndarray) -> np.ndarray:
    """
    Normalize confidence scores to 0-1 range
    
    Args:
        scores: Raw confidence scores
        
    Returns:
        Normalized scores
    """
    if not scores:
        return []
    
    # 转换为列表处理（如果输入是numpy数组）
    scores_list = list(scores)
    
    # 处理空列表
    if len(scores_list) == 0:
        return []
    
    # 找到最大值进行标准化
    max_score = max(scores_list)
    
    if max_score == 0:
        return [0.0] * len(scores_list)
    
    normalized = [score / max_score for score in scores_list]
    return normalized