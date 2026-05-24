"""
NumberedSequence class for antibody/nanobody sequence representation.

Stores ANARCII-numbered sequences with FR/CDR decomposition.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from panig.scheme import (
    CDRBoundaries,
    SchemeDefinition,
    classify_position,
    get_scheme,
)


@dataclass
class NumberedPosition:
    """A single numbered position in an antibody sequence."""
    position: int          # Numbered position (e.g., IMGT 1-128)
    residue: str           # Amino acid at this position
    region: str            # FR1, CDR1, FR2, CDR2, FR3, CDR3, or FR4


@dataclass
class NumberedSequence:
    """
    A numbered antibody or nanobody sequence with FR/CDR decomposition.

    Attributes:
        name: Sequence identifier
        sequence: Original raw sequence
        chain_type: 'heavy' (VH), 'nanobody' (VHH), or 'light' (VL)
        scheme: Numbering scheme used (e.g., 'imgt', 'kabat')
        positions: List of NumberedPosition objects
        species: Detected species (from ANARCII germline assignment)
        v_gene: V-gene assignment
        j_gene: J-gene assignment
        v_gene_identity: V-gene percent identity
    """
    name: str
    sequence: str
    chain_type: str  # 'heavy', 'nanobody', 'light'
    scheme: str
    positions: List[NumberedPosition] = field(default_factory=list)
    species: Optional[str] = None
    v_gene: Optional[str] = None
    j_gene: Optional[str] = None
    v_gene_identity: Optional[float] = None

    @property
    def scheme_def(self) -> SchemeDefinition:
        """Get the SchemeDefinition for this sequence's numbering scheme."""
        return get_scheme(self.scheme)

    @property
    def cdr1(self) -> str:
        """Get CDR1 sequence."""
        return self._get_region("CDR1")

    @property
    def cdr2(self) -> str:
        """Get CDR2 sequence."""
        return self._get_region("CDR2")

    @property
    def cdr3(self) -> str:
        """Get CDR3 sequence."""
        return self._get_region("CDR3")

    @property
    def fr1(self) -> str:
        """Get FR1 sequence."""
        return self._get_region("FR1")

    @property
    def fr2(self) -> str:
        """Get FR2 sequence."""
        return self._get_region("FR2")

    @property
    def fr3(self) -> str:
        """Get FR3 sequence."""
        return self._get_region("FR3")

    @property
    def fr4(self) -> str:
        """Get FR4 sequence."""
        return self._get_region("FR4")

    def _get_region(self, region: str) -> str:
        """Get concatenated residues for a specific region."""
        return "".join(
            p.residue for p in self.positions if p.region == region
        )

    def get_region_positions(self, region: str) -> List[NumberedPosition]:
        """Get all positions belonging to a specific region."""
        return [p for p in self.positions if p.region == region]

    def get_framework_positions(self) -> List[NumberedPosition]:
        """Get all framework positions (FR1+FR2+FR3+FR4)."""
        return [
            p for p in self.positions
            if p.region.startswith("FR")
        ]

    def get_cdr_positions(self) -> List[NumberedPosition]:
        """Get all CDR positions (CDR1+CDR2+CDR3)."""
        return [
            p for p in self.positions
            if p.region.startswith("CDR")
        ]

    def get_full_sequence(self) -> str:
        """Get the full numbered sequence (all positions)."""
        return "".join(p.residue for p in self.positions)

    def get_framework_sequence(self) -> str:
        """Get concatenated framework sequence (FR1+FR2+FR3+FR4)."""
        return "".join(
            p.residue for p in self.positions
            if p.region.startswith("FR")
        )

    def to_fasta(self) -> str:
        """Export as FASTA format."""
        return f">{self.name}\n{self.sequence}"

    def summary(self) -> Dict[str, str]:
        """Get a summary dictionary of the sequence regions."""
        return {
            "name": self.name,
            "chain_type": self.chain_type,
            "scheme": self.scheme,
            "species": self.species or "unknown",
            "v_gene": self.v_gene or "unknown",
            "FR1": self.fr1,
            "CDR1": self.cdr1,
            "FR2": self.fr2,
            "CDR2": self.cdr2,
            "FR3": self.fr3,
            "CDR3": self.cdr3,
            "FR4": self.fr4,
        }

    def __len__(self) -> int:
        return len(self.positions)

    def __str__(self) -> str:
        return (
            f"NumberedSequence(name={self.name}, "
            f"chain={self.chain_type}, "
            f"scheme={self.scheme}, "
            f"length={len(self)})"
        )
