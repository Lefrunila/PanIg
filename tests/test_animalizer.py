"""Tests for the xenotypizer module."""

import pytest
from panig.species_profiles import SpeciesProfile
from panig.xenotypizer import Xenotypizer, XenotypizationResult


class TestSpeciesProfile:
    """Test SpeciesProfile class."""

    def test_create_profile(self):
        """Test creating a species profile."""
        profile = SpeciesProfile(
            species="dog",
            chain_type="VH",
            profile={
                1: {"E": 0.8, "D": 0.2},
                2: {"V": 0.9, "L": 0.1},
            }
        )

        assert profile.species == "dog"
        assert profile.chain_type == "VH"
        assert len(profile.profile) == 2

    def test_get_frequency(self):
        """Test getting amino acid frequency."""
        profile = SpeciesProfile(
            species="dog",
            chain_type="VH",
            profile={
                1: {"E": 0.8, "D": 0.2},
            }
        )

        assert profile.get_frequency(1, "E") == 0.8
        assert profile.get_frequency(1, "D") == 0.2
        assert profile.get_frequency(1, "A") == 0.0
        assert profile.get_frequency(999, "E") == 0.0

    def test_get_consensus(self):
        """Test getting consensus amino acid."""
        profile = SpeciesProfile(
            species="dog",
            chain_type="VH",
            profile={
                1: {"E": 0.8, "D": 0.2},
                2: {"V": 0.4, "L": 0.6},
            }
        )

        assert profile.get_consensus(1) == "E"
        assert profile.get_consensus(2) == "L"
        assert profile.get_consensus(999) is None

    def test_is_native(self):
        """Test native residue detection."""
        profile = SpeciesProfile(
            species="dog",
            chain_type="VH",
            profile={
                1: {"E": 0.8, "D": 0.2},
                2: {"V": 0.05, "L": 0.95},
            }
        )

        # E at position 1 is native (frequency 0.8 >= 0.1)
        assert profile.is_native(1, "E", threshold=0.1) is True
        # D at position 1 is native (frequency 0.2 >= 0.1)
        assert profile.is_native(1, "D", threshold=0.1) is True
        # A at position 1 is not native (frequency 0.0 < 0.1)
        assert profile.is_native(1, "A", threshold=0.1) is False

    def test_get_substitution(self):
        """Test substitution suggestion."""
        profile = SpeciesProfile(
            species="dog",
            chain_type="VH",
            profile={
                1: {"E": 0.8, "D": 0.2},
                2: {"V": 0.05, "L": 0.95},
            }
        )

        # No substitution needed for native residue
        assert profile.get_substitution(1, "E", threshold=0.1) is None
        # Substitution needed for non-native residue
        assert profile.get_substitution(2, "V", threshold=0.1) == "L"

    def test_save_load(self, tmp_path):
        """Test saving and loading profiles."""
        profile = SpeciesProfile(
            species="dog",
            chain_type="VH",
            profile={
                1: {"E": 0.8, "D": 0.2},
            }
        )

        # Save
        save_path = tmp_path / "test_profile.json"
        profile.save(str(save_path))

        # Load
        loaded = SpeciesProfile.load(str(save_path))
        assert loaded.species == "dog"
        assert loaded.chain_type == "VH"
        assert loaded.get_frequency(1, "E") == 0.8

    def test_from_sequences(self):
        """Test building profile from sequences."""
        from panig.sequence import NumberedPosition, NumberedSequence

        # Create mock sequences
        positions1 = [
            NumberedPosition(1, "E", "FR1"),
            NumberedPosition(2, "V", "FR1"),
        ]
        positions2 = [
            NumberedPosition(1, "E", "FR1"),
            NumberedPosition(2, "L", "FR1"),
        ]

        seq1 = NumberedSequence(
            name="seq1",
            sequence="EV",
            chain_type="heavy",
            scheme="imgt",
            positions=positions1,
        )
        seq2 = NumberedSequence(
            name="seq2",
            sequence="EL",
            chain_type="heavy",
            scheme="imgt",
            positions=positions2,
        )

        profile = SpeciesProfile.from_sequences(
            species="test",
            chain_type="VH",
            sequences=[seq1, seq2],
        )

        assert profile.get_frequency(1, "E") == 1.0
        assert profile.get_frequency(2, "V") == 0.5
        assert profile.get_frequency(2, "L") == 0.5

    def test_get_compatible_substitution(self):
        """Test compatible substitution with constrained allowed set."""
        profile = SpeciesProfile(
            species="dog",
            chain_type="VH",
            profile={
                1: {"E": 0.75, "D": 0.15, "K": 0.05, "A": 0.03, "V": 0.02},
            }
        )
        # A at position 1 is not native; only D and E are allowed (ionic negative)
        result = profile.get_compatible_substitution(1, "A", {"D", "E"}, threshold=0.1)
        assert result == "E"  # E has highest freq among allowed

    def test_get_compatible_substitution_already_native(self):
        """Test that no substitution is returned when current AA is already native."""
        profile = SpeciesProfile(
            species="dog",
            chain_type="VH",
            profile={
                1: {"E": 0.75, "D": 0.15, "K": 0.05, "A": 0.03},
            }
        )
        # E at position 1 is native (0.75 >= 0.1)
        result = profile.get_compatible_substitution(1, "E", {"D", "E"}, threshold=0.1)
        assert result is None

    def test_get_compatible_substitution_no_compatible_candidates(self):
        """Test when no compatible AA meets the threshold."""
        profile = SpeciesProfile(
            species="dog",
            chain_type="VH",
            profile={
                1: {"E": 0.75, "D": 0.15, "K": 0.05, "A": 0.03},
            }
        )
        # A is not native, but only {"K", "R"} are allowed, and K is only 0.05 (< 0.1)
        result = profile.get_compatible_substitution(1, "A", {"K", "R"}, threshold=0.1)
        assert result is None


class TestXenotypizer:
    """Test Xenotypizer class."""

    def test_xenotypizer_initialization(self):
        """Test xenotypizer initialization."""
        xenotypizer = Xenotypizer(threshold=0.15)
        assert xenotypizer.threshold == 0.15

    def test_xenotypization_result(self):
        """Test XenotypizationResult dataclass."""
        from panig.sequence import NumberedPosition, NumberedSequence

        positions = [
            NumberedPosition(1, "E", "FR1"),
            NumberedPosition(2, "V", "FR1"),
            NumberedPosition(30, "S", "CDR1"),
        ]

        original = NumberedSequence(
            name="test",
            sequence="EV...S",
            chain_type="nanobody",
            scheme="imgt",
            positions=positions,
        )

        result = XenotypizationResult(
            original=original,
            modified_sequence="DL...S",
            modified_name="test_xenotypized_dog",
            target_species="dog",
            operation="xenotypize",
        )

        assert result.modified_name == "test_xenotypized_dog"
        assert result.target_species == "dog"
        assert result.operation == "xenotypize"

    def test_humanization_result(self):
        """Test XenotypizationResult for humanization."""
        from panig.sequence import NumberedPosition, NumberedSequence

        positions = [
            NumberedPosition(1, "E", "FR1"),
            NumberedPosition(2, "V", "FR1"),
            NumberedPosition(30, "S", "CDR1"),
        ]

        original = NumberedSequence(
            name="test",
            sequence="EV...S",
            chain_type="nanobody",
            scheme="imgt",
            positions=positions,
        )

        result = XenotypizationResult(
            original=original,
            modified_sequence="DL...S",
            modified_name="test_humanized",
            target_species="human",
            operation="humanize",
        )

        assert result.modified_name == "test_humanized"
        assert result.target_species == "human"
        assert result.operation == "humanize"

    def test_backward_compatibility(self):
        """Test backward compatibility aliases."""
        from panig.sequence import NumberedPosition, NumberedSequence

        positions = [
            NumberedPosition(1, "E", "FR1"),
            NumberedPosition(2, "V", "FR1"),
        ]

        original = NumberedSequence(
            name="test",
            sequence="EV",
            chain_type="nanobody",
            scheme="imgt",
            positions=positions,
        )

        result = XenotypizationResult(
            original=original,
            modified_sequence="DL",
            modified_name="test_xenotypized",
            target_species="dog",
            operation="xenotypize",
        )

        # Test backward compatibility aliases
        assert result.xenotypized_sequence == "DL"
        assert result.xenotypized_name == "test_xenotypized"
        assert result.animalized_sequence == "DL"  # backward compat
        assert result.animalized_name == "test_xenotypized"  # backward compat
