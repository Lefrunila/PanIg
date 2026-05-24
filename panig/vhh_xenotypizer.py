"""VHH (nanobody) xenotypization module."""
from __future__ import annotations

import csv
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

from panig.scheme import get_scheme, classify_position
from panig.sequence import NumberedPosition, NumberedSequence
from panig.species_profiles import SpeciesProfile

logger = logging.getLogger(__name__)


# VHH-specific IMGT positions that must be preserved
# These positions have evolved for single-domain stability
# and must NOT be changed to conventional VH residues
VHH_LOCKED_POSITIONS = {
    37,   # FR2-CDR1 junction
    44,   # Former VH-VL interface
    45,   # Former VH-VL interface
    47,   # VHH hallmark - tryptophan solubility anchor
    83,   # FR3 core packing
    84,   # Charged pocket (Glu46-Lys82-Lys104 triad)
    103,  # CDR3 base aromatic
    108,  # FR4 stability
}


@dataclass
class VHHSubstitutionRecord:
    """Record of a single amino acid substitution in VHH."""
    position: int
    region: str
    original_aa: str
    new_aa: str
    original_frequency: float
    new_frequency: float
    excluded: bool = False
    exclusion_reason: Optional[str] = None
    interaction_constrained: bool = False


@dataclass
class VHHXenotypizationResult:
    """Result of a VHH xenotypization or humanization operation."""
    original: NumberedSequence
    modified_sequence: str
    modified_name: str
    target_species: str
    operation: str  # 'xenotypize' or 'humanize'
    substitutions: List[VHHSubstitutionRecord] = field(default_factory=list)
    excluded_positions: Set[int] = field(default_factory=set)
    interaction_map: Dict[int, str] = field(default_factory=dict)
    locked_positions: Set[int] = field(default_factory=set)

    @property
    def total_substitutions(self) -> int:
        """Number of positions that were changed."""
        return len([s for s in self.substitutions if not s.excluded and s.original_aa != s.new_aa])

    @property
    def total_excluded(self) -> int:
        """Number of positions excluded from substitution."""
        return len([s for s in self.substitutions if s.excluded])

    def summary(self) -> Dict[str, any]:
        """Get a summary of the result."""
        return {
            "original_name": self.original.name,
            "target_species": self.target_species,
            "operation": self.operation,
            "total_positions": len(self.substitutions),
            "substitutions": self.total_substitutions,
            "excluded": self.total_excluded,
            "locked_vhh": len(self.locked_positions),
            "chain_type": self.original.chain_type,
            "scheme": self.original.scheme,
        }

    def to_csv(self, path: str):
        """Export substitution report to CSV."""
        with open(path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "Position", "Region", "Original_AA", "New_AA",
                "Original_Freq", "New_Freq", "Excluded", "Reason",
                "Interaction_Constrained",
            ])
            for sub in self.substitutions:
                writer.writerow([
                    sub.position,
                    sub.region,
                    sub.original_aa,
                    sub.new_aa,
                    f"{sub.original_frequency:.4f}",
                    f"{sub.new_frequency:.4f}",
                    sub.excluded,
                    sub.exclusion_reason or "",
                    sub.interaction_constrained,
                ])


class VHHXenotypizer:
    """
    VHH (nanobody) xenotypization engine.

    Mirrors Llamanade's approach but adapted for xenotypization:
    - Uses target species frequency profiles
    - Preserves VHH-specific framework positions
    - Preserves CDR residues
    - Excludes structurally interacting residues
    """

    def __init__(
        self,
        threshold: float = 0.1,
        scheme: str = "imgt",
        use_protinter: bool = True,
    ):
        """
        Initialize the VHH xenotypizer.

        Args:
            threshold: Frequency threshold for determining native residues
            scheme: Numbering scheme to use
            use_protinter: Whether to use Protinter for interaction detection
        """
        self.threshold = threshold
        self.scheme = scheme
        self.scheme_def = get_scheme(scheme)
        self.use_protinter = use_protinter

    def xenotypize(
        self,
        numbered_seq: NumberedSequence,
        target_species: str,
        species_profile: Optional[SpeciesProfile] = None,
        pdb_file: Optional[str] = None,
        interaction_map: Optional[Dict[int, str]] = None,
    ) -> VHHXenotypizationResult:
        """
        Xenotypize a VHH nanobody to a target species.

        Args:
            numbered_seq: Numbered VHH sequence
            target_species: Target species name
            species_profile: Pre-loaded species frequency profile
            pdb_file: PDB file for interaction detection (optional)
            interaction_map: Pre-computed interaction map (pos -> type).
                If provided, skips interaction detection. For testing.

        Returns:
            VHHXenotypizationResult with the xenotypized sequence
        """
        from panig.interactions import get_compatible_residues

        # Load species profile if not provided
        if species_profile is None:
            species_profile = self._load_species_profile(
                target_species, numbered_seq.chain_type
            )

        # Get structural interactions as interaction map
        if interaction_map is None:
            if self.use_protinter:
                interaction_map = self._detect_interactions(
                    numbered_seq, pdb_file=pdb_file
                )
            else:
                interaction_map = {}

        # Process each position
        substitutions = []
        new_residues = {}

        for pos in numbered_seq.positions:
            # Skip CDR positions - never modify
            if pos.region.startswith("CDR"):
                substitutions.append(VHHSubstitutionRecord(
                    position=pos.position,
                    region=pos.region,
                    original_aa=pos.residue,
                    new_aa=pos.residue,
                    original_frequency=species_profile.get_frequency(
                        pos.position, pos.residue
                    ),
                    new_frequency=species_profile.get_frequency(
                        pos.position, pos.residue
                    ),
                    excluded=True,
                    exclusion_reason="cdr",
                ))
                continue

            # Check if position is VHH-specific (locked)
            if pos.position in VHH_LOCKED_POSITIONS:
                substitutions.append(VHHSubstitutionRecord(
                    position=pos.position,
                    region=pos.region,
                    original_aa=pos.residue,
                    new_aa=pos.residue,
                    original_frequency=species_profile.get_frequency(
                        pos.position, pos.residue
                    ),
                    new_frequency=species_profile.get_frequency(
                        pos.position, pos.residue
                    ),
                    excluded=True,
                    exclusion_reason="vhh_specific",
                ))
                continue

            current_freq = species_profile.get_frequency(
                pos.position, pos.residue
            )

            # Check if position is involved in structural interactions
            if pos.position in interaction_map:
                int_type = interaction_map[pos.position]
                compatible = get_compatible_residues(pos.residue, int_type)
                suggestion = species_profile.get_compatible_substitution(
                    pos.position, pos.residue, compatible, self.threshold
                )

                if suggestion is not None:
                    new_freq = species_profile.get_frequency(
                        pos.position, suggestion
                    )
                    substitutions.append(VHHSubstitutionRecord(
                        position=pos.position,
                        region=pos.region,
                        original_aa=pos.residue,
                        new_aa=suggestion,
                        original_frequency=current_freq,
                        new_frequency=new_freq,
                        excluded=False,
                        interaction_constrained=True,
                    ))
                    new_residues[pos.position] = suggestion
                else:
                    substitutions.append(VHHSubstitutionRecord(
                        position=pos.position,
                        region=pos.region,
                        original_aa=pos.residue,
                        new_aa=pos.residue,
                        original_frequency=current_freq,
                        new_frequency=current_freq,
                        excluded=True,
                        exclusion_reason="structural_interaction",
                    ))
                continue

            # Check if residue needs substitution
            if current_freq >= self.threshold:
                # Already native to target species
                substitutions.append(VHHSubstitutionRecord(
                    position=pos.position,
                    region=pos.region,
                    original_aa=pos.residue,
                    new_aa=pos.residue,
                    original_frequency=current_freq,
                    new_frequency=current_freq,
                    excluded=False,
                ))
            else:
                # Needs substitution - use consensus
                consensus = species_profile.get_consensus(pos.position)
                if consensus:
                    new_freq = species_profile.get_frequency(
                        pos.position, consensus
                    )
                    substitutions.append(VHHSubstitutionRecord(
                        position=pos.position,
                        region=pos.region,
                        original_aa=pos.residue,
                        new_aa=consensus,
                        original_frequency=current_freq,
                        new_frequency=new_freq,
                        excluded=False,
                    ))
                    new_residues[pos.position] = consensus
                else:
                    # No consensus available - keep original
                    substitutions.append(VHHSubstitutionRecord(
                        position=pos.position,
                        region=pos.region,
                        original_aa=pos.residue,
                        new_aa=pos.residue,
                        original_frequency=current_freq,
                        new_frequency=current_freq,
                        excluded=False,
                    ))

        # Build the modified sequence
        modified_seq = self._build_modified_sequence(
            numbered_seq, new_residues
        )

        # Create result
        result = VHHXenotypizationResult(
            original=numbered_seq,
            modified_sequence=modified_seq,
            modified_name=f"{numbered_seq.name}_xenotypized_{target_species}",
            target_species=target_species,
            operation="xenotypize",
            substitutions=substitutions,
            excluded_positions=set(interaction_map.keys()),
            interaction_map=interaction_map,
            locked_positions=VHH_LOCKED_POSITIONS,
        )

        logger.info(
            f"Xenotypized {numbered_seq.name} to {target_species}: "
            f"{result.total_substitutions} substitutions, "
            f"{result.total_excluded} excluded, "
            f"{len(result.locked_positions)} VHH-locked"
        )

        return result

    # Backward compatibility alias
    animalize = xenotypize

    def humanize(
        self,
        numbered_seq: NumberedSequence,
        species_profile: Optional[SpeciesProfile] = None,
        pdb_file: Optional[str] = None,
    ) -> VHHXenotypizationResult:
        """
        Humanize a VHH nanobody.

        Args:
            numbered_seq: Numbered VHH sequence
            species_profile: Pre-loaded human species profile
            pdb_file: PDB file for interaction detection (optional)

        Returns:
            VHHXenotypizationResult with the humanized sequence
        """
        return self.xenotypize(
            numbered_seq=numbered_seq,
            target_species="human",
            species_profile=species_profile,
            pdb_file=pdb_file,
        )

    def _load_species_profile(
        self, species: str, chain_type: str
    ) -> SpeciesProfile:
        """Load species profile from built-in profiles or cache."""

        chain_map = {"heavy": "VH", "nanobody": "VHH", "light": "VL"}
        file_suffix = chain_map.get(chain_type, chain_type)

        # Build list of suffixes to try (VHH -> VH fallback for nanobodies)
        suffixes_to_try = [file_suffix]
        if file_suffix == "VHH":
            suffixes_to_try.append("VH")

        profile_dir = Path(__file__).parent.parent / "profiles"
        cache_dir = Path.home() / ".panig" / "cache" / "profiles"

        for suffix in suffixes_to_try:
            # Try local profiles directory
            profile_file = profile_dir / f"{species}_{suffix}.json"
            if profile_file.exists():
                return SpeciesProfile.load(str(profile_file))

            # Try Google Drive cache
            cache_file = cache_dir / f"{species}_{suffix}.json"
            if cache_file.exists():
                return SpeciesProfile.load(str(cache_file))

        raise FileNotFoundError(
            f"Species profile not found for {species} ({chain_type}). "
            f"Run 'panig download --species {species}' first."
        )

    def _detect_interactions(
        self, numbered_seq: NumberedSequence, pdb_file: Optional[str] = None
    ) -> Dict[int, str]:
        """
        Detect structural interactions using ImmuneBuilder + Protinter.

        If a PDB file is provided, uses it directly. Otherwise predicts
        structure via ImmuneBuilder. Falls back to sequence-based heuristics
        if structure prediction fails.

        Args:
            numbered_seq: Numbered VHH sequence
            pdb_file: Optional path to existing PDB file

        Returns:
            Dict mapping position -> interaction type
        """
        from panig.interactions import InteractionDetector
        from panig.structure import StructurePredictor

        detector = InteractionDetector()

        # If no PDB file provided, predict structure via ImmuneBuilder
        if pdb_file is None:
            try:
                predictor = StructurePredictor()
                seq = "".join(p.residue for p in numbered_seq.positions)
                pdb_file = predictor.predict(
                    seq, chain_type="nanobody", name=numbered_seq.name
                )
            except Exception as e:
                logger.warning(
                    f"Structure prediction failed, using sequence heuristics: {e}"
                )
                interactions = detector.detect_from_sequence(
                    sequence=None, numbered_positions=numbered_seq.positions
                )
                return detector.get_interaction_map(interactions)

        # Use PDB-based detection
        try:
            interactions = detector.detect_interactions(pdb_path=pdb_file)
        except Exception as e:
            logger.warning(
                f"PDB interaction detection failed, using sequence heuristics: {e}"
            )
            interactions = detector.detect_from_sequence(
                sequence=None, numbered_positions=numbered_seq.positions
            )

        return detector.get_interaction_map(interactions)

    def _build_modified_sequence(
        self,
        numbered_seq: NumberedSequence,
        substitutions: Dict[int, str],
    ) -> str:
        """Build the modified sequence with substitutions applied."""
        result = []
        for pos in numbered_seq.positions:
            if pos.position in substitutions:
                result.append(substitutions[pos.position])
            else:
                result.append(pos.residue)
        return "".join(result)


# Backward compatibility aliases
VHHAnimalizationResult = VHHXenotypizationResult
VHHAnimalizer = VHHXenotypizer
