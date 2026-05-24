"""
Species-specific antibody/nanobody frequency profiles.

Manages pre-computed amino acid frequency profiles for different species,
used to determine which framework residues need adaptation.
"""

import json
import logging
from pathlib import Path
from typing import Dict, Optional, Set

import numpy as np

logger = logging.getLogger(__name__)


class SpeciesProfile:
    """
    Species-specific amino acid frequency profile for antibody frameworks.

    Each profile is a dictionary mapping position -> {amino_acid: frequency}.
    Used to determine which residues are "native" to a target species.
    """

    # Standard amino acids
    AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

    def __init__(
        self,
        species: str,
        chain_type: str = "VH",
        profile: Optional[Dict[int, Dict[str, float]]] = None,
    ):
        """
        Initialize a species profile.

        Args:
            species: Species name (e.g., 'human', 'dog', 'cat')
            chain_type: Chain type ('VH' for heavy, 'VHH' for nanobody)
            profile: Pre-loaded profile data {position: {aa: frequency}}
        """
        self.species = species.lower()
        self.chain_type = chain_type
        self.profile: Dict[int, Dict[str, float]] = profile or {}

    def get_frequency(self, position: int, amino_acid: str) -> float:
        """
        Get the frequency of an amino acid at a specific position.

        Args:
            position: Numbered position (e.g., IMGT position)
            amino_acid: Single-letter amino acid code

        Returns:
            Frequency (0.0 to 1.0), or 0.0 if position/AA not in profile
        """
        pos_profile = self.profile.get(position, {})
        return pos_profile.get(amino_acid.upper(), 0.0)

    def get_consensus(self, position: int) -> Optional[str]:
        """
        Get the most common amino acid at a position.

        Args:
            position: Numbered position

        Returns:
            Most frequent amino acid, or None if position not in profile
        """
        pos_profile = self.profile.get(position)
        if not pos_profile:
            return None
        return max(pos_profile, key=pos_profile.get)

    def get_consensus_frequency(self, position: int) -> float:
        """
        Get the frequency of the consensus amino acid at a position.

        Args:
            position: Numbered position

        Returns:
            Frequency of the most common amino acid (0.0 to 1.0)
        """
        pos_profile = self.profile.get(position)
        if not pos_profile:
            return 0.0
        return max(pos_profile.values())

    def is_native(
        self,
        position: int,
        amino_acid: str,
        threshold: float = 0.1,
    ) -> bool:
        """
        Check if an amino acid is 'native' to this species at a position.

        A residue is considered native if its frequency in the species
        profile is above the threshold.

        Args:
            position: Numbered position
            amino_acid: Single-letter amino acid code
            threshold: Minimum frequency to consider native (default: 0.1)

        Returns:
            True if the amino acid is common in this species at this position
        """
        return self.get_frequency(position, amino_acid) >= threshold

    def get_substitution(
        self,
        position: int,
        current_aa: str,
        threshold: float = 0.1,
    ) -> Optional[str]:
        """
        Get a species-appropriate substitution for a position.

        Only suggests substitution if the current residue is below threshold.

        Args:
            position: Numbered position
            current_aa: Current amino acid at this position
            threshold: Minimum frequency to consider native

        Returns:
            Suggested substitution amino acid, or None if no change needed
        """
        if self.is_native(position, current_aa, threshold):
            return None
        return self.get_consensus(position)

    def get_compatible_substitution(
        self,
        position: int,
        current_aa: str,
        allowed_aas: Set[str],
        threshold: float = 0.1,
    ) -> Optional[str]:
        """
        Get a substitution from a chemically constrained set of amino acids.

        Returns the highest-frequency amino acid from allowed_aas that
        differs from current_aa and meets the frequency threshold.
        Returns None if current_aa is already above threshold among the
        allowed set, or if no allowed AA meets the threshold.

        Args:
            position: Numbered position
            current_aa: Current amino acid at this position
            allowed_aas: Set of amino acids allowed by chemical compatibility
            threshold: Minimum frequency to consider

        Returns:
            Best compatible substitution, or None
        """
        pos_profile = self.profile.get(position, {})

        # If the current residue is already above threshold, no substitution needed
        if pos_profile.get(current_aa, 0.0) >= threshold:
            return None

        # Find the highest-frequency AA from the allowed set (excluding current)
        candidates = {
            aa: freq for aa, freq in pos_profile.items()
            if aa in allowed_aas and aa != current_aa and freq >= threshold
        }

        if not candidates:
            return None

        return max(candidates, key=candidates.get)

    def save(self, path: str):
        """Save profile to a JSON file."""
        data = {
            "species": self.species,
            "chain_type": self.chain_type,
            "profile": {str(k): v for k, v in self.profile.items()},
        }
        with open(path, "w") as f:
            json.dump(data, f, indent=2)

    @classmethod
    def load(cls, path: str) -> "SpeciesProfile":
        """
        Load a profile from a JSON file.

        Args:
            path: Path to the JSON profile file

        Returns:
            SpeciesProfile instance
        """
        with open(path, "r") as f:
            data = json.load(f)

        profile = {int(k): v for k, v in data["profile"].items()}

        return cls(
            species=data["species"],
            chain_type=data.get("chain_type", "VH"),
            profile=profile,
        )

    @classmethod
    def from_sequences(
        cls,
        species: str,
        chain_type: str,
        sequences: list,
        scheme: str = "imgt",
    ) -> "SpeciesProfile":
        """
        Build a frequency profile from a list of numbered sequences.

        Args:
            species: Species name
            chain_type: Chain type
            sequences: List of NumberedSequence objects
            scheme: Numbering scheme

        Returns:
            SpeciesProfile built from the sequences
        """
        from collections import defaultdict

        # Count amino acids at each position
        position_counts: Dict[int, Dict[str, int]] = defaultdict(
            lambda: defaultdict(int)
        )

        for seq in sequences:
            for pos in seq.positions:
                position_counts[pos.position][pos.residue] += 1

        # Convert counts to frequencies
        profile = {}
        for position, counts in position_counts.items():
            total = sum(counts.values())
            if total > 0:
                profile[position] = {
                    aa: count / total
                    for aa, count in counts.items()
                }

        return cls(
            species=species,
            chain_type=chain_type,
            profile=profile,
        )

    def __repr__(self) -> str:
        return (
            f"SpeciesProfile(species={self.species}, "
            f"chain_type={self.chain_type}, "
            f"positions={len(self.profile)})"
        )
