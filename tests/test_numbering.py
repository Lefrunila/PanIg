"""Tests for the numbering module."""

import pytest
from panig.numbering import Numberer
from panig.sequence import NumberedSequence
from panig.scheme import get_scheme, classify_position


class TestScheme:
    """Test numbering scheme definitions."""

    def test_get_scheme_imgt(self):
        """Test getting IMGT scheme."""
        scheme = get_scheme("imgt")
        assert scheme.name == "IMGT"
        assert scheme.cdr_boundaries.cdr1 == (27, 38)

    def test_get_scheme_kabat(self):
        """Test getting Kabat scheme."""
        scheme = get_scheme("kabat")
        assert scheme.name == "Kabat"

    def test_get_scheme_case_insensitive(self):
        """Test case-insensitive scheme lookup."""
        assert get_scheme("IMGT").name == "IMGT"
        assert get_scheme("Imgt").name == "IMGT"

    def test_get_scheme_invalid(self):
        """Test invalid scheme raises error."""
        with pytest.raises(ValueError, match="Unknown numbering scheme"):
            get_scheme("invalid")

    def test_classify_position_cdr1(self):
        """Test CDR1 classification."""
        scheme = get_scheme("imgt")
        assert classify_position(30, scheme) == "CDR1"
        assert classify_position(27, scheme) == "CDR1"
        assert classify_position(38, scheme) == "CDR1"

    def test_classify_position_fr1(self):
        """Test FR1 classification."""
        scheme = get_scheme("imgt")
        assert classify_position(1, scheme) == "FR1"
        assert classify_position(26, scheme) == "FR1"

    def test_classify_position_cdr3(self):
        """Test CDR3 classification."""
        scheme = get_scheme("imgt")
        assert classify_position(105, scheme) == "CDR3"
        assert classify_position(117, scheme) == "CDR3"


class TestNumberedSequence:
    """Test NumberedSequence class."""

    def test_create_sequence(self):
        """Test creating a numbered sequence."""
        from panig.sequence import NumberedPosition

        positions = [
            NumberedPosition(1, "E", "FR1"),
            NumberedPosition(2, "V", "FR1"),
            NumberedPosition(30, "S", "CDR1"),
            NumberedPosition(105, "C", "CDR3"),
        ]

        seq = NumberedSequence(
            name="test",
            sequence="EV...S...C",
            chain_type="nanobody",
            scheme="imgt",
            positions=positions,
        )

        assert seq.name == "test"
        assert seq.chain_type == "nanobody"
        assert seq.fr1 == "EV"
        assert seq.cdr1 == "S"
        assert seq.cdr3 == "C"

    def test_get_framework_positions(self):
        """Test getting framework positions."""
        from panig.sequence import NumberedPosition

        positions = [
            NumberedPosition(1, "E", "FR1"),
            NumberedPosition(2, "V", "FR1"),
            NumberedPosition(30, "S", "CDR1"),
            NumberedPosition(105, "C", "CDR3"),
            NumberedPosition(120, "G", "FR4"),
        ]

        seq = NumberedSequence(
            name="test",
            sequence="EV...S...C...G",
            chain_type="nanobody",
            scheme="imgt",
            positions=positions,
        )

        fr_positions = seq.get_framework_positions()
        assert len(fr_positions) == 3
        assert all(p.region.startswith("FR") for p in fr_positions)


class TestNumberer:
    """Test the Numberer class."""

    def test_parse_fasta(self, tmp_path):
        """Test FASTA parsing."""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">seq1\nACDEFG\n>seq2\nHIKLMN\n")

        sequences = Numberer._parse_fasta(str(fasta_file))
        assert len(sequences) == 2
        assert sequences["seq1"] == "ACDEFG"
        assert sequences["seq2"] == "HIKLMN"

    def test_parse_fasta_multiline(self, tmp_path):
        """Test FASTA parsing with multiline sequences."""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">seq1\nACDEFG\nHIKLMN\n")

        sequences = Numberer._parse_fasta(str(fasta_file))
        assert sequences["seq1"] == "ACDEFGHIKLMN"
