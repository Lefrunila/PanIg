"""Tests for the scorer module."""

import shutil
from pathlib import Path

import pytest

from panig.scorer import Scorer


# Check if blastp is available
_BLASTP_AVAILABLE = shutil.which("blastp") is not None


class TestScorer:
    """Test Scorer class."""

    def test_scorer_initialization(self):
        """Test scorer initialization."""
        scorer = Scorer(blastdb_path="/tmp/test_db")
        assert scorer.blastdb_path == "/tmp/test_db"
        assert scorer.evalue == 0.001

    def test_calculate_t20(self):
        """Test T20 score calculation."""
        scorer = Scorer()

        # Test with hits
        hits = [
            {"bitscore": 100, "pident": 90.0},
            {"bitscore": 90, "pident": 85.0},
            {"bitscore": 80, "pident": 80.0},
        ]

        score = scorer._calculate_t20(hits)
        assert score == pytest.approx(85.0)  # (90 + 85 + 80) / 3

    def test_calculate_t20_empty(self):
        """Test T20 score with no hits."""
        scorer = Scorer()
        assert scorer._calculate_t20([]) == 0.0

    def test_calculate_t20_sorts_by_bitscore(self):
        """Test that T20 sorts hits by bitscore descending before averaging."""
        scorer = Scorer()

        hits = [
            {"bitscore": 50, "pident": 70.0},
            {"bitscore": 200, "pident": 99.0},
            {"bitscore": 100, "pident": 85.0},
        ]

        # Average of top 3 by bitscore: 99.0, 85.0, 70.0 = 84.666...
        score = scorer._calculate_t20(hits)
        assert score == pytest.approx((99.0 + 85.0 + 70.0) / 3)

    def test_calculate_t20_single_hit(self):
        """Test T20 with a single hit."""
        scorer = Scorer()
        hits = [{"bitscore": 100, "pident": 92.5}]
        assert scorer._calculate_t20(hits) == pytest.approx(92.5)

    def test_calculate_t20_max_20_hits(self):
        """Test that T20 uses at most 20 hits."""
        scorer = Scorer()

        # Create 30 hits with varying pident
        hits = [
            {"bitscore": 300 - i, "pident": 95.0 - i}
            for i in range(30)
        ]

        score = scorer._calculate_t20(hits)

        # Should average only top 20
        expected = sum(95.0 - i for i in range(20)) / 20
        assert score == pytest.approx(expected)

    def test_parse_blast_output(self):
        """Test BLAST output parsing."""
        scorer = Scorer()

        output = """# BLASTP 2.12.0+
# Query: test
# Database: test_db
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
test\tsubject1\t90.00\t100\t10\t0\t1\t100\t1\t100\t1e-50\t200
test\tsubject2\t85.00\t95\t14\t0\t1\t95\t1\t95\t1e-45\t180
# 2 hits found
"""

        hits = scorer._parse_blast_output(output)
        assert len(hits) == 2
        assert hits[0]["pident"] == 90.0
        assert hits[0]["bitscore"] == 200.0
        assert hits[0]["subject_id"] == "subject1"
        assert hits[0]["length"] == 100
        assert hits[0]["mismatch"] == 10
        assert hits[0]["gapopen"] == 0
        assert hits[0]["qstart"] == 1
        assert hits[0]["qend"] == 100
        assert hits[0]["sstart"] == 1
        assert hits[0]["send"] == 100
        assert hits[0]["evalue"] == pytest.approx(1e-50)
        assert hits[1]["pident"] == 85.0

    def test_parse_blast_output_empty(self):
        """Test parsing empty BLAST output."""
        scorer = Scorer()
        assert scorer._parse_blast_output("") == []
        assert scorer._parse_blast_output("# No hits found\n") == []

    def test_parse_blast_output_malformed_lines(self):
        """Test that malformed lines are skipped gracefully."""
        scorer = Scorer()

        output = """# BLASTP 2.12.0+
test\tsubject1\t90.00\t100\t10\t0\t1\t100\t1\t100\t1e-50\t200
bad_line_with_few_fields
test\tsubject2\t85.00\t95\t14\t0\t1\t95\t1\t95\t1e-45\t180
# 2 hits found
"""

        hits = scorer._parse_blast_output(output)
        assert len(hits) == 2
        assert hits[0]["subject_id"] == "subject1"
        assert hits[1]["subject_id"] == "subject2"

    def test_score_sequence_no_db_raises(self):
        """Test that score_sequence raises ValueError when no DB is configured."""
        scorer = Scorer()
        with pytest.raises(ValueError, match="No BLAST database configured"):
            scorer.score_sequence("ACDEFG")


@pytest.mark.skipif(
    not _BLASTP_AVAILABLE,
    reason="blastp not installed",
)
class TestScorerIntegration:
    """Integration tests using a real BLAST database."""

    @pytest.fixture(autouse=True)
    def _setup_scorer(self):
        blastdb_path = Path(__file__).parent.parent / "blastdb" / "dog_VH_blastdb" / "dog_VH"
        self.scorer = Scorer(blastdb_path=str(blastdb_path))
        self.vhh_sequence = (
            "QVQLVESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAITWSGGNTYY"
            "ADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAADRGYYGSGYWGQGTQVTVSS"
        )

    def test_score_sequence_returns_positive(self):
        """Test that scoring a VHH sequence returns a positive T20 score."""
        score = self.scorer.score_sequence(self.vhh_sequence)
        assert isinstance(score, float)
        assert score > 0.0, f"Expected positive T20 score, got {score}"

    def test_score_sequence_range(self):
        """Test that T20 score is in valid range (0-100)."""
        score = self.scorer.score_sequence(self.vhh_sequence)
        assert 0.0 <= score <= 100.0
