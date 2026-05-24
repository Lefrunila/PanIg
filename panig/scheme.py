"""
Numbering schemes and CDR/FR boundary definitions for antibodies and nanobodies.

Supports IMGT, Kabat, Chothia, Martin, and AHo numbering schemes.
Based on Llamanade's Scheme.py with extensions for multiple chain types.
"""

from dataclasses import dataclass
from typing import Dict, List, Tuple


@dataclass(frozen=True)
class CDRBoundaries:
    """CDR loop boundaries for a specific numbering scheme."""
    cdr1: Tuple[int, int]  # (start, end) positions
    cdr2: Tuple[int, int]
    cdr3: Tuple[int, int]


@dataclass(frozen=True)
class SchemeDefinition:
    """Complete definition of a numbering scheme."""
    name: str
    cdr_boundaries: CDRBoundaries
    fr_boundaries: Dict[str, Tuple[int, int]]  # FR1, FR2, FR3, FR4


# IMGT numbering scheme (international standard)
IMGT_CDR = CDRBoundaries(
    cdr1=(27, 38),
    cdr2=(56, 65),
    cdr3=(105, 117),
)

IMGT_FR = {
    "FR1": (1, 26),
    "FR2": (39, 55),
    "FR3": (66, 104),
    "FR4": (118, 128),
}

IMGT = SchemeDefinition(
    name="IMGT",
    cdr_boundaries=IMGT_CDR,
    fr_boundaries=IMGT_FR,
)


# Kabat numbering scheme
KABAT_CDR = CDRBoundaries(
    cdr1=(31, 35),
    cdr2=(50, 65),
    cdr3=(95, 102),
)

KABAT_FR = {
    "FR1": (1, 30),
    "FR2": (36, 49),
    "FR3": (66, 94),
    "FR4": (103, 113),
}

KABAT = SchemeDefinition(
    name="Kabat",
    cdr_boundaries=KABAT_CDR,
    fr_boundaries=KABAT_FR,
)


# Chothia numbering scheme
CHOTHIA_CDR = CDRBoundaries(
    cdr1=(26, 32),
    cdr2=(52, 56),
    cdr3=(95, 102),
)

CHOTHIA_FR = {
    "FR1": (1, 25),
    "FR2": (33, 51),
    "FR3": (57, 94),
    "FR4": (103, 113),
}

CHOTHIA = SchemeDefinition(
    name="Chothia",
    cdr_boundaries=CHOTHIA_CDR,
    fr_boundaries=CHOTHIA_FR,
)


# Martin (extended Chothia) numbering scheme
MARTIN_CDR = CDRBoundaries(
    cdr1=(24, 34),
    cdr2=(50, 56),
    cdr3=(89, 101),
)

MARTIN_FR = {
    "FR1": (1, 23),
    "FR2": (35, 49),
    "FR3": (57, 88),
    "FR4": (102, 113),
}

MARTIN = SchemeDefinition(
    name="Martin",
    cdr_boundaries=MARTIN_CDR,
    fr_boundaries=MARTIN_FR,
)


# AHo numbering scheme
AHO_CDR = CDRBoundaries(
    cdr1=(25, 40),
    cdr2=(56, 69),
    cdr3=(105, 117),
)

AHO_FR = {
    "FR1": (1, 24),
    "FR2": (41, 55),
    "FR3": (70, 104),
    "FR4": (118, 128),
}

AHO = SchemeDefinition(
    name="AHo",
    cdr_boundaries=AHO_CDR,
    fr_boundaries=AHO_FR,
)


# Registry of all schemes
SCHEMES: Dict[str, SchemeDefinition] = {
    "imgt": IMGT,
    "kabat": KABAT,
    "chothia": CHOTHIA,
    "martin": MARTIN,
    "aho": AHO,
}


def get_scheme(name: str) -> SchemeDefinition:
    """Get a numbering scheme by name (case-insensitive)."""
    scheme = SCHEMES.get(name.lower())
    if scheme is None:
        raise ValueError(
            f"Unknown numbering scheme: {name}. "
            f"Available: {', '.join(SCHEMES.keys())}"
        )
    return scheme


def classify_position(
    position: int,
    scheme: SchemeDefinition,
) -> str:
    """
    Classify a numbered position as FR1, CDR1, FR2, CDR2, FR3, CDR3, or FR4.

    Args:
        position: The numbered position (1-indexed).
        scheme: The numbering scheme to use.

    Returns:
        Region name (e.g., "FR1", "CDR1", etc.)
    """
    # Check CDRs first
    for i, (start, end) in enumerate([
        scheme.cdr_boundaries.cdr1,
        scheme.cdr_boundaries.cdr2,
        scheme.cdr_boundaries.cdr3,
    ], 1):
        if start <= position <= end:
            return f"CDR{i}"

    # Check FRs
    for fr_name, (start, end) in scheme.fr_boundaries.items():
        if start <= position <= end:
            return fr_name

    # Position outside defined regions
    return "Unknown"


def get_framework_positions(scheme: SchemeDefinition) -> List[int]:
    """Get all framework positions for a given scheme."""
    positions = []
    for fr_name, (start, end) in scheme.fr_boundaries.items():
        positions.extend(range(start, end + 1))
    return sorted(positions)


def get_cdr_positions(scheme: SchemeDefinition) -> List[int]:
    """Get all CDR positions for a given scheme."""
    positions = []
    for i, (start, end) in enumerate([
        scheme.cdr_boundaries.cdr1,
        scheme.cdr_boundaries.cdr2,
        scheme.cdr_boundaries.cdr3,
    ], 1):
        positions.extend(range(start, end + 1))
    return sorted(positions)
