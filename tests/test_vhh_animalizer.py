"""Tests for the VHH animalizer module."""

import pytest

from panig.scheme import classify_position, get_scheme
from panig.sequence import NumberedPosition, NumberedSequence
from panig.species_profiles import SpeciesProfile
from panig.vhh_xenotypizer import VHHXenotypizer, VHHXenotypizationResult, VHH_LOCKED_POSITIONS


# VHH sequence used for testing
VHH_SEQUENCE = (
    "QVQLVESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAITWSGGNTYY"
    "ADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAADRGYYGSGYWGQGTQVTVSS"
)


def _build_numbered_vhh(sequence: str = VHH_SEQUENCE) -> NumberedSequence:
    """
    Build a NumberedSequence for a VHH using IMGT numbering.

    Maps each residue to its IMGT position (1..128) and classifies
    the region (FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4).
    """
    scheme = get_scheme("imgt")

    imgt_positions = list(range(1, len(sequence) + 1))

    positions = []
    for i, res in enumerate(sequence):
        pos_num = imgt_positions[i]
        region = classify_position(pos_num, scheme)
        positions.append(NumberedPosition(
            position=pos_num,
            residue=res,
            region=region,
        ))

    return NumberedSequence(
        name="test_vhh",
        sequence=sequence,
        chain_type="nanobody",
        scheme="imgt",
        positions=positions,
    )


def _build_dog_profile(sequence: str = VHH_SEQUENCE) -> SpeciesProfile:
    """
    Build a realistic mock dog species profile.

    For each position, the actual VHH residue gets high frequency (0.75)
    so it won't be substituted. A second "dog-preferred" residue gets
    moderate frequency (0.15) and the rest get small baseline values.
    This ensures:
    - Most framework residues are already "native" (above threshold)
    - A few positions may get substituted (simulating real animalization)
    - Sequence identity stays >80%
    """
    scheme = get_scheme("imgt")
    amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
    profile_data = {}

    for i, res in enumerate(sequence):
        pos_num = i + 1
        freqs = {}

        # Baseline for all AAs
        for aa in amino_acids:
            freqs[aa] = 0.005

        # The actual VHH residue gets high frequency
        freqs[res] = 0.75

        # A different "dog consensus" residue gets moderate frequency
        dog_preferred = amino_acids[(i * 7 + 3) % 20]
        if dog_preferred != res:
            freqs[dog_preferred] = 0.15

        profile_data[pos_num] = freqs

    return SpeciesProfile(
        species="dog",
        chain_type="VHH",
        profile=profile_data,
    )


class TestVHHXenotypizer:
    """Test VHHXenotypizer class."""

    def test_instantiation(self):
        """Test VHHXenotypizer instantiation."""
        xenotypizer = VHHXenotypizer()
        assert xenotypizer.threshold == 0.1
        assert xenotypizer.scheme == "imgt"
        assert xenotypizer.use_protinter is True

    def test_instantiation_custom_params(self):
        """Test VHHXenotypizer with custom parameters."""
        xenotypizer = VHHXenotypizer(
            threshold=0.15,
            scheme="kabat",
            use_protinter=False,
        )
        assert xenotypizer.threshold == 0.15
        assert xenotypizer.scheme == "kabat"
        assert xenotypizer.use_protinter is False

    def test_vhh_locked_positions_defined(self):
        """Test that VHH-locked positions are defined correctly."""
        expected = {37, 44, 45, 47, 83, 84, 103, 108}
        assert VHH_LOCKED_POSITIONS == expected

    def test_xenotypize_returns_result(self):
        """Test that xenotypize returns a VHHXenotypizationResult."""
        xenotypizer = VHHXenotypizer(threshold=0.1, use_protinter=False)
        numbered_seq = _build_numbered_vhh()
        profile = _build_dog_profile()

        result = xenotypizer.xenotypize(
            numbered_seq=numbered_seq,
            target_species="dog",
            species_profile=profile,
        )

        assert isinstance(result, VHHXenotypizationResult)
        assert result.target_species == "dog"
        assert result.operation == "xenotypize"
        assert result.original is numbered_seq

    def test_xenotypized_sequence_length_matches_original(self):
        """Test that the xenotypized sequence has the same length as the original."""
        xenotypizer = VHHXenotypizer(threshold=0.1, use_protinter=False)
        numbered_seq = _build_numbered_vhh()
        profile = _build_dog_profile()

        result = xenotypizer.xenotypize(
            numbered_seq=numbered_seq,
            target_species="dog",
            species_profile=profile,
        )

        assert len(result.modified_sequence) == len(VHH_SEQUENCE)

    def test_cdr_positions_preserved(self):
        """Test that CDR residues are never changed during xenotypization."""
        xenotypizer = VHHXenotypizer(threshold=0.0, use_protinter=False)
        numbered_seq = _build_numbered_vhh()
        profile = _build_dog_profile()

        result = xenotypizer.xenotypize(
            numbered_seq=numbered_seq,
            target_species="dog",
            species_profile=profile,
        )

        # Check that all CDR positions have excluded=True with reason "cdr"
        cdr_subs = [
            s for s in result.substitutions
            if s.region.startswith("CDR")
        ]

        for sub in cdr_subs:
            assert sub.excluded is True, (
                f"CDR position {sub.position} should be excluded"
            )
            assert sub.exclusion_reason == "cdr", (
                f"CDR position {sub.position} should have reason 'cdr', "
                f"got '{sub.exclusion_reason}'"
            )
            # The new AA should be the same as the original
            assert sub.new_aa == sub.original_aa, (
                f"CDR position {sub.position}: original={sub.original_aa}, "
                f"new={sub.new_aa}"
            )

    def test_vhh_locked_positions_preserved(self):
        """Test that VHH-locked positions are preserved during xenotypization."""
        xenotypizer = VHHXenotypizer(threshold=0.0, use_protinter=False)
        numbered_seq = _build_numbered_vhh()
        profile = _build_dog_profile()

        result = xenotypizer.xenotypize(
            numbered_seq=numbered_seq,
            target_species="dog",
            species_profile=profile,
        )

        # Check positions 44, 45, 47, 83, 84, 103
        for pos_num in [44, 45, 47, 83, 84, 103]:
            matching = [
                s for s in result.substitutions if s.position == pos_num
            ]
            if matching:
                sub = matching[0]
                assert sub.excluded is True, (
                    f"VHH-locked position {pos_num} should be excluded"
                )
                assert sub.exclusion_reason == "vhh_specific", (
                    f"Position {pos_num} should have reason 'vhh_specific', "
                    f"got '{sub.exclusion_reason}'"
                )
                assert sub.new_aa == sub.original_aa, (
                    f"VHH-locked position {pos_num}: original={sub.original_aa}, "
                    f"new={sub.new_aa}"
                )

    def test_sequence_identity_reasonable(self):
        """Test that xenotypized sequence has reasonable identity to original (>80%)."""
        xenotypizer = VHHXenotypizer(threshold=0.1, use_protinter=False)
        numbered_seq = _build_numbered_vhh()
        profile = _build_dog_profile()

        result = xenotypizer.xenotypize(
            numbered_seq=numbered_seq,
            target_species="dog",
            species_profile=profile,
        )

        # Calculate sequence identity
        original = VHH_SEQUENCE
        modified = result.modified_sequence
        matches = sum(
            1 for a, b in zip(original, modified) if a == b
        )
        identity = matches / len(original) * 100

        assert identity > 80.0, (
            f"Sequence identity {identity:.1f}% is below 80% threshold"
        )

    def test_result_summary(self):
        """Test that result summary returns expected fields."""
        xenotypizer = VHHXenotypizer(threshold=0.1, use_protinter=False)
        numbered_seq = _build_numbered_vhh()
        profile = _build_dog_profile()

        result = xenotypizer.xenotypize(
            numbered_seq=numbered_seq,
            target_species="dog",
            species_profile=profile,
        )

        summary = result.summary()

        assert summary["target_species"] == "dog"
        assert summary["operation"] == "xenotypize"
        assert summary["chain_type"] == "nanobody"
        assert summary["scheme"] == "imgt"
        assert "total_positions" in summary
        assert "substitutions" in summary
        assert "excluded" in summary
        assert "locked_vhh" in summary

    def test_result_locked_positions_set(self):
        """Test that the result records the VHH-locked positions."""
        xenotypizer = VHHXenotypizer(threshold=0.1, use_protinter=False)
        numbered_seq = _build_numbered_vhh()
        profile = _build_dog_profile()

        result = xenotypizer.xenotypize(
            numbered_seq=numbered_seq,
            target_species="dog",
            species_profile=profile,
        )

        assert result.locked_positions == VHH_LOCKED_POSITIONS

    def test_humanize_delegates_to_xenotypize(self):
        """Test that humanize() delegates to xenotypize with species='human'."""
        xenotypizer = VHHXenotypizer(threshold=0.1, use_protinter=False)
        numbered_seq = _build_numbered_vhh()

        # Build a human profile
        human_profile = SpeciesProfile(
            species="human",
            chain_type="VHH",
            profile=_build_dog_profile().profile,  # Reuse structure
        )

        result = xenotypizer.humanize(
            numbered_seq=numbered_seq,
            species_profile=human_profile,
        )

        assert result.target_species == "human"
        assert result.operation == "xenotypize"  # humanize delegates to xenotypize
