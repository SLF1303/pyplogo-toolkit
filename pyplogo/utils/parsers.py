"""
Utility functions for parsing and handling protein data
"""

def parse_fasta(fasta_string: str) -> str:
    """
    Parse FASTA string and extract sequence
    
    Args:
        fasta_string: FASTA format string
        
    Returns:
        Protein sequence as string
        
    Raises:
        ValueError: If input is empty, contains only whitespace, 
                   or has no sequence data after headers
    """
    if not isinstance(fasta_string, str):
        raise TypeError(f"Expected string input, got {type(fasta_string).__name__}")
    
    stripped_input = fasta_string.strip()
    if not stripped_input:
        raise ValueError("FASTA input cannot be empty or contain only whitespace")
    
    lines = stripped_input.split('\n')
    sequence_lines = []
    
    for line in lines:
        stripped_line = line.strip()
        if stripped_line and not stripped_line.startswith('>'):
            sequence_lines.append(stripped_line)
    
    sequence = ''.join(sequence_lines)
    
    if not sequence:
        raise ValueError("No sequence data found in FASTA input")
    
    return sequence

def validate_sequence(sequence: str) -> bool:
    """
    Validate amino acid sequence
    
    Args:
        sequence: Amino acid sequence to validate
        
    Returns:
        True if valid, False otherwise
    """
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    return all(aa in valid_aa for aa in sequence)

def calculate_sequence_stats(sequence: str) -> dict:
    """
    Calculate basic statistics for amino acid sequence
    
    Args:
        sequence: Amino acid sequence
        
    Returns:
        Dictionary with sequence statistics
    """
    from collections import Counter
    from ..data.aa_properties import AA_PROPERTIES
    
    counter = Counter(sequence)
    total = len(sequence)
    
    stats = {
        'length': total,
        'aa_composition': {aa: count/total for aa, count in counter.items()},
        'type_composition': {},
        'molecular_weight': 0.0
    }
    
    # Calculate type composition and molecular weight
    type_counter = {}
    for aa, count in counter.items():
        aa_type = AA_PROPERTIES.get(aa, {}).get('type', 'unknown')
        type_counter[aa_type] = type_counter.get(aa_type, 0) + count
        stats['molecular_weight'] += count * AA_PROPERTIES.get(aa, {}).get('weight', 0)
    
    stats['type_composition'] = {t: count/total for t, count in type_counter.items()}
    
    return stats