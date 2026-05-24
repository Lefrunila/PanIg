"""
Utility functions for PanIg.

I/O helpers, FASTA parsing, and common operations.
"""

import csv
import logging
from pathlib import Path
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)


def parse_fasta(fasta_path: str) -> Dict[str, str]:
    """
    Parse a FASTA file into a dictionary of {name: sequence}.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Dictionary mapping sequence names to sequences
    """
    sequences = {}
    current_name = None
    current_seq = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = "".join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    if current_name is not None:
        sequences[current_name] = "".join(current_seq)

    return sequences


def write_fasta(sequences: Dict[str, str], output_path: str):
    """
    Write sequences to a FASTA file.

    Args:
        sequences: Dictionary of {name: sequence}
        output_path: Path to output FASTA file
    """
    with open(output_path, "w") as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n{seq}\n")


def validate_sequence(sequence: str) -> Tuple[bool, str]:
    """
    Validate an amino acid sequence.

    Args:
        sequence: Amino acid sequence string

    Returns:
        Tuple of (is_valid, error_message)
    """
    valid_chars = set("ACDEFGHIKLMNPQRSTVWY")

    if not sequence:
        return False, "Empty sequence"

    invalid_chars = set(sequence.upper()) - valid_chars
    if invalid_chars:
        return False, f"Invalid characters: {invalid_chars}"

    if len(sequence) < 50:
        return False, f"Sequence too short ({len(sequence)} aa, minimum 50)"

    if len(sequence) > 200:
        return False, f"Sequence too long ({len(sequence)} aa, maximum 200)"

    return True, ""


def format_sequence(sequence: str, width: int = 60) -> str:
    """
    Format a sequence with line breaks for readability.

    Args:
        sequence: Amino acid sequence
        width: Line width

    Returns:
        Formatted sequence string
    """
    return "\n".join(
        sequence[i:i + width] for i in range(0, len(sequence), width)
    )


def create_output_directory(path: str) -> Path:
    """
    Create an output directory if it doesn't exist.

    Args:
        path: Directory path

    Returns:
        Path object
    """
    output_dir = Path(path)
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def load_csv(path: str) -> List[Dict[str, str]]:
    """
    Load a CSV file into a list of dictionaries.

    Args:
        path: Path to CSV file

    Returns:
        List of row dictionaries
    """
    rows = []
    with open(path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def save_csv(data: List[Dict[str, str]], path: str, fieldnames: List[str] = None):
    """
    Save a list of dictionaries to CSV.

    Args:
        data: List of row dictionaries
        path: Path to output CSV file
        fieldnames: Column names (auto-detected if None)
    """
    if not data:
        return

    if fieldnames is None:
        fieldnames = list(data[0].keys())

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)


def get_sequence_stats(sequence: str) -> Dict[str, any]:
    """
    Get basic statistics for a sequence.

    Args:
        sequence: Amino acid sequence

    Returns:
        Dictionary of statistics
    """
    if not sequence:
        return {}

    length = len(sequence)
    aa_counts = {}
    for aa in sequence:
        aa_counts[aa] = aa_counts.get(aa, 0) + 1

    # Calculate molecular weight (approximate)
    mw_table = {
        'A': 89.09, 'C': 121.16, 'D': 133.10, 'E': 147.13,
        'F': 165.19, 'G': 75.03, 'H': 155.16, 'I': 131.17,
        'K': 146.19, 'L': 131.17, 'M': 149.21, 'N': 132.12,
        'P': 115.13, 'Q': 146.15, 'R': 174.20, 'S': 105.09,
        'T': 119.12, 'V': 117.15, 'W': 204.23, 'Y': 181.19,
    }
    mw = sum(mw_table.get(aa, 0) for aa in sequence)

    return {
        "length": length,
        "molecular_weight": round(mw, 2),
        "aa_composition": aa_counts,
    }
