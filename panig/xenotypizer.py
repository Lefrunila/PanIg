"""
Core xenotypization logic.

Adapts antibody/nanobody framework regions to match target species preferences
while preserving CDR sequences and structurally important residues.

Terminology:
- Xenotypize: General term for adapting to a non-native species (inverse of humanize)
- Species-specific: caninize (dog), felinize (cat), equinize (horse), bovinize (cattle)
"""

import csv
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

from panig.numbering import Numberer
from panig.scheme import get_scheme
from panig.sequence import NumberedSequence
from panig.species_profiles import SpeciesProfile

logger = logging.getLogger(__name__)


@dataclass
class SubstitutionRecord:
    """Record of a single amino acid substitution."""
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
class XenotypizationResult:
    """Result of a xenotypization or humanization operation."""
    original: NumberedSequence
    modified_sequence: str
    modified_name: str
    target_species: str
    operation: str  # 'xenotypize' or 'humanize'
    substitutions: List[SubstitutionRecord] = field(default_factory=list)
    excluded_positions: Set[int] = field(default_factory=set)
    interaction_map: Dict[int, str] = field(default_factory=dict)

    @property
    def xenotypized_sequence(self) -> str:
        """Alias for modified_sequence."""
        return self.modified_sequence

    @property
    def xenotypized_name(self) -> str:
        """Alias for modified_name."""
        return self.modified_name

    # Backward compatibility aliases
    @property
    def animalized_sequence(self) -> str:
        """Alias for modified_sequence (backward compatibility)."""
        return self.modified_sequence

    @property
    def animalized_name(self) -> str:
        """Alias for modified_name (backward compatibility)."""
        return self.modified_name

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


class Xenotypizer:
    """
    Antibody/nanobody xenotypization and humanization engine.

    Supports two operations:
    - Xenotypize: Adapt sequences to match a target species
    - Humanize: Adapt sequences to match human preferences

    Preserves:
    - All CDR residues (antigen-binding regions)
    - Structurally important residues (intra-chain interactions)
    """

    def __init__(
        self,
        numberer: Optional[Numberer] = None,
        threshold: float = 0.1,
        scheme: str = "imgt",
        use_synthetic: bool = False,
    ):
        """
        Initialize the xenotypizer.

        Args:
            numberer: Numberer instance for sequence numbering
            threshold: Frequency threshold for determining native residues
            scheme: Numbering scheme to use
            use_synthetic: If True, prefer synthetic profiles over germline-only
                         profiles for species with sparse data (e.g. cat/goat VL).
                         Experimental — synthetic profiles may not improve results.
        """
        self.numberer = numberer or Numberer(scheme=scheme)
        self.threshold = threshold
        self.scheme = scheme
        self.use_synthetic = use_synthetic

    def xenotypize(
        self,
        sequence: str,
        target_species: str,
        name: str = "input",
        chain_type: Optional[str] = None,
        interacting_positions: Optional[Union[Set[int], Dict[int, str]]] = None,
        species_profile: Optional[SpeciesProfile] = None,
    ) -> XenotypizationResult:
        """
        Xenotypize an antibody/nanobody sequence to a target species.

        Args:
            sequence: Input amino acid sequence
            target_species: Target species name (e.g., 'dog', 'cat')
            name: Sequence identifier
            chain_type: Chain type ('heavy', 'nanobody', 'light')
            interacting_positions: Positions involved in structural
                interactions. Accepts Set[int] (all treated as ionic)
                or Dict[int, str] mapping position -> interaction type.
            species_profile: Pre-loaded species profile. If None, will try
                           to load from built-in profiles.

        Returns:
            XenotypizationResult with the xenotypized sequence and metadata
        """
        # Number the input sequence
        numbered = self.numberer.number_sequence(sequence, name, chain_type)

        # Load species profile if not provided
        if species_profile is None:
            species_profile = self._load_species_profile(
                target_species, numbered.chain_type
            )

        # Get interaction map if not provided
        if interacting_positions is None:
            interaction_map = self._detect_interactions(numbered)
        elif isinstance(interacting_positions, dict):
            interaction_map = interacting_positions
        else:
            # Legacy Set[int] input — treat all as ionic (safest)
            interaction_map = {pos: "ionic" for pos in interacting_positions}

        # Perform xenotypization
        substitutions, new_residues = self._substitute_residues(
            numbered, species_profile, interaction_map
        )

        # Build the modified sequence
        modified_seq = self._build_modified_sequence(numbered, new_residues)

        # Create result
        result = XenotypizationResult(
            original=numbered,
            modified_sequence=modified_seq,
            modified_name=f"{name}_xenotypized_{target_species}",
            target_species=target_species,
            operation="xenotypize",
            substitutions=substitutions,
            excluded_positions=set(interaction_map.keys()),
            interaction_map=interaction_map,
        )

        logger.info(
            f"Xenotypized {name} to {target_species}: "
            f"{result.total_substitutions} substitutions, "
            f"{result.total_excluded} excluded"
        )

        return result

    def humanize(
        self,
        sequence: str,
        name: str = "input",
        chain_type: Optional[str] = None,
        interacting_positions: Optional[Union[Set[int], Dict[int, str]]] = None,
        species_profile: Optional[SpeciesProfile] = None,
    ) -> XenotypizationResult:
        """
        Humanize an antibody/nanobody sequence.

        This is equivalent to xenotypizing to human species.
        Used to reduce immunogenicity in human therapeutic applications.

        Args:
            sequence: Input amino acid sequence
            name: Sequence identifier
            chain_type: Chain type ('heavy', 'nanobody', 'light')
            interacting_positions: Positions involved in structural
                interactions. Accepts Set[int] or Dict[int, str].
            species_profile: Pre-loaded human species profile. If None, will try
                           to load from built-in profiles.

        Returns:
            XenotypizationResult with the humanized sequence and metadata
        """
        # Number the input sequence
        numbered = self.numberer.number_sequence(sequence, name, chain_type)

        # Load human profile if not provided
        if species_profile is None:
            species_profile = self._load_species_profile(
                "human", numbered.chain_type
            )

        # Get interaction map if not provided
        if interacting_positions is None:
            interaction_map = self._detect_interactions(numbered)
        elif isinstance(interacting_positions, dict):
            interaction_map = interacting_positions
        else:
            interaction_map = {pos: "ionic" for pos in interacting_positions}

        # Perform humanization
        substitutions, new_residues = self._substitute_residues(
            numbered, species_profile, interaction_map
        )

        # Build the modified sequence
        modified_seq = self._build_modified_sequence(numbered, new_residues)

        # Create result
        result = XenotypizationResult(
            original=numbered,
            modified_sequence=modified_seq,
            modified_name=f"{name}_humanized",
            target_species="human",
            operation="humanize",
            substitutions=substitutions,
            excluded_positions=set(interaction_map.keys()),
            interaction_map=interaction_map,
        )

        logger.info(
            f"Humanized {name}: "
            f"{result.total_substitutions} substitutions, "
            f"{result.total_excluded} excluded"
        )

        return result

    def _substitute_residues(
        self,
        numbered: NumberedSequence,
        species_profile: SpeciesProfile,
        interaction_map: Dict[int, str],
    ) -> Tuple[List[SubstitutionRecord], Dict[int, str]]:
        """
        Perform residue substitutions based on species profile.

        Args:
            numbered: Numbered input sequence
            species_profile: Target species profile
            interaction_map: Positions -> interaction type mapping

        Returns:
            Tuple of (substitutions list, new residues dict)
        """
        from panig.interactions import get_compatible_residues

        substitutions = []
        new_residues = {}

        for pos in numbered.positions:
            # Skip CDR positions - never modify
            if pos.region.startswith("CDR"):
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
                    substitutions.append(SubstitutionRecord(
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
                    substitutions.append(SubstitutionRecord(
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
            suggestion = species_profile.get_substitution(
                pos.position, pos.residue, self.threshold
            )

            if suggestion is not None:
                # Substitution needed
                new_freq = species_profile.get_frequency(
                    pos.position, suggestion
                )
                substitutions.append(SubstitutionRecord(
                    position=pos.position,
                    region=pos.region,
                    original_aa=pos.residue,
                    new_aa=suggestion,
                    original_frequency=current_freq,
                    new_frequency=new_freq,
                    excluded=False,
                ))
                new_residues[pos.position] = suggestion
            else:
                # No substitution needed
                substitutions.append(SubstitutionRecord(
                    position=pos.position,
                    region=pos.region,
                    original_aa=pos.residue,
                    new_aa=pos.residue,
                    original_frequency=current_freq,
                    new_frequency=current_freq,
                    excluded=False,
                ))

        return substitutions, new_residues

    def xenotypize_batch(
        self,
        sequences: Dict[str, str],
        target_species: str,
        chain_type: Optional[str] = None,
        species_profile: Optional[SpeciesProfile] = None,
    ) -> List[XenotypizationResult]:
        """
        Xenotypize multiple sequences in batch.

        Args:
            sequences: Dictionary of {name: sequence} pairs
            target_species: Target species name
            chain_type: Chain type for all sequences
            species_profile: Pre-loaded species profile

        Returns:
            List of XenotypizationResult objects
        """
        results = []
        for name, seq in sequences.items():
            try:
                result = self.xenotypize(
                    seq, target_species, name, chain_type, species_profile=species_profile
                )
                results.append(result)
            except Exception as e:
                logger.warning(f"Failed to xenotypize {name}: {e}")
                continue
        return results

    def humanize_batch(
        self,
        sequences: Dict[str, str],
        chain_type: Optional[str] = None,
        species_profile: Optional[SpeciesProfile] = None,
    ) -> List[XenotypizationResult]:
        """
        Humanize multiple sequences in batch.

        Args:
            sequences: Dictionary of {name: sequence} pairs
            chain_type: Chain type for all sequences
            species_profile: Pre-loaded human species profile

        Returns:
            List of XenotypizationResult objects
        """
        results = []
        for name, seq in sequences.items():
            try:
                result = self.humanize(
                    seq, name, chain_type, species_profile=species_profile
                )
                results.append(result)
            except Exception as e:
                logger.warning(f"Failed to humanize {name}: {e}")
                continue
        return results

    def _load_species_profile(
        self, species: str, chain_type: str
    ) -> SpeciesProfile:
        """
        Load a species profile from built-in profiles or Google Drive cache.

        For nanobodies (VHH), falls back to VH profile if no VHH-specific
        profile exists, since nanobody framework regions are similar to VH.

        Args:
            species: Species name
            chain_type: Chain type ('heavy', 'nanobody', 'light')

        Returns:
            SpeciesProfile instance

        Raises:
            FileNotFoundError: If profile not found
        """
        # Map chain_type to file suffix
        chain_map = {"heavy": "VH", "nanobody": "VHH", "light": "VL"}
        file_suffix = chain_map.get(chain_type, chain_type)

        # Build list of suffixes to try (VHH -> VH fallback for nanobodies)
        suffixes_to_try = [file_suffix]
        if file_suffix == "VHH":
            suffixes_to_try.append("VH")  # Fallback to VH profile

        profile_dir = Path(__file__).parent.parent / "profiles"
        cache_dir = Path.home() / ".panig" / "cache" / "profiles"

        for suffix in suffixes_to_try:
            # Try synthetic profile first if requested
            if self.use_synthetic:
                synth_file = profile_dir / f"{species}_{suffix}_synthetic.json"
                if synth_file.exists():
                    logger.info(
                        f"Using synthetic profile for {species} ({suffix})"
                    )
                    return SpeciesProfile.load(str(synth_file))

            # Try local profiles directory
            profile_file = profile_dir / f"{species}_{suffix}.json"
            if profile_file.exists():
                if suffix != file_suffix:
                    logger.info(
                        f"Using {suffix} profile for {species} ({chain_type})"
                    )
                return SpeciesProfile.load(str(profile_file))

            # Try Google Drive cache
            cache_file = cache_dir / f"{species}_{suffix}.json"
            if cache_file.exists():
                if suffix != file_suffix:
                    logger.info(
                        f"Using {suffix} profile for {species} ({chain_type})"
                    )
                return SpeciesProfile.load(str(cache_file))

        # Try downloading from Google Drive
        try:
            self._download_profile(species, chain_type)
            cache_file = cache_dir / f"{species}_{file_suffix}.json"
            if cache_file.exists():
                return SpeciesProfile.load(str(cache_file))
        except Exception as e:
            logger.warning(f"Could not download profile: {e}")

        raise FileNotFoundError(
            f"Species profile not found for {species} ({chain_type}). "
            f"Tried: {', '.join(f'{species}_{s}.json' for s in suffixes_to_try)}\n"
            f"Run 'panig download --species {species}' to fetch from Google Drive."
        )

    def _download_profile(self, species: str, chain_type: str):
        """Download a species profile from Google Drive."""
        import subprocess

        chain_map = {"heavy": "VH", "nanobody": "VHH", "light": "VL"}
        file_suffix = chain_map.get(chain_type, chain_type)

        cache_dir = Path.home() / ".panig" / "cache" / "profiles"
        cache_dir.mkdir(parents=True, exist_ok=True)

        remote_path = f"gdrive:PanIg_databases/profiles/{species}_{file_suffix}.json"
        local_path = str(cache_dir / f"{species}_{file_suffix}.json")

        subprocess.run(
            ["rclone", "copy", remote_path, str(cache_dir)],
            check=True,
            capture_output=True,
        )

    def _detect_interactions(
        self, numbered: NumberedSequence
    ) -> Dict[int, str]:
        """
        Detect structural interactions using ImmuneBuilder + Protinter.

        Predicts 3D structure via ImmuneBuilder, then uses InteractionDetector
        to identify intra-chain interactions (ionic, cation-pi, pi-pi).
        Falls back to sequence-based heuristics if structure prediction fails.

        Args:
            numbered: Numbered sequence

        Returns:
            Dict mapping position -> interaction type
        """
        from panig.interactions import InteractionDetector
        from panig.structure import StructurePredictor

        detector = InteractionDetector()

        # Try structure-based detection via ImmuneBuilder
        try:
            predictor = StructurePredictor()
            chain_type = "nanobody" if numbered.chain_type == "nanobody" else "antibody"
            seq = "".join(p.residue for p in numbered.positions)
            pdb_path = predictor.predict(
                seq, chain_type=chain_type, name=numbered.name
            )
            interactions = detector.detect_interactions(pdb_path=pdb_path)
        except Exception as e:
            logger.warning(
                f"Structure prediction failed, using sequence heuristics: {e}"
            )
            interactions = detector.detect_from_sequence(
                sequence=None, numbered_positions=numbered.positions
            )

        return detector.get_interaction_map(interactions)

    def _build_modified_sequence(
        self,
        numbered: NumberedSequence,
        substitutions: Dict[int, str],
    ) -> str:
        """
        Build the modified sequence with substitutions applied.

        Args:
            numbered: Original numbered sequence
            substitutions: Dict of {position: new_residue}

        Returns:
            Modified amino acid sequence
        """
        result = []
        for pos in numbered.positions:
            if pos.position in substitutions:
                result.append(substitutions[pos.position])
            else:
                result.append(pos.residue)
        return "".join(result)


# Backward compatibility aliases
AnimalizationResult = XenotypizationResult
Animalizer = Xenotypizer
