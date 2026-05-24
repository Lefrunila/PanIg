"""
Species-specific scoring for xenotypized antibodies/nanobodies.

Implements T20-like scoring that measures how well a sequence
matches the target species' antibody repertoire.
"""

import csv
import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from panig.sequence import NumberedSequence
from panig.species_profiles import SpeciesProfile

logger = logging.getLogger(__name__)


class Scorer:
    """
    Species-specific antibody/nanobody scoring.

    Computes a T20-like score by BLASTing framework regions against
    a target species antibody database. Higher scores indicate
    better species match.
    """

    def __init__(
        self,
        blastdb_path: Optional[str] = None,
        evalue: float = 0.001,
        num_threads: int = 4,
    ):
        """
        Initialize the scorer.

        Args:
            blastdb_path: Path to species-specific BLAST database.
                         Can be a directory (auto-detects .pin/.psq prefix)
                         or a file prefix (e.g., 'dog_VH_blastdb/dog_VH').
            evalue: E-value threshold for BLAST
            num_threads: Number of threads for BLAST
        """
        self.blastdb_path = self._resolve_blastdb_path(blastdb_path)
        self.evalue = evalue
        self.num_threads = num_threads

    @staticmethod
    def _resolve_blastdb_path(path: Optional[str]) -> Optional[str]:
        """Resolve BLAST database path to a file prefix.

        BLAST+ needs the file prefix (e.g., 'db/dog_VH'), not the
        directory ('db/dog_VH_blastdb/'). This auto-detects the prefix.
        """
        if path is None:
            return None

        p = Path(path)

        # Check if it's a directory
        if p.is_dir():
            pins = list(p.glob("*.pin"))
            if pins:
                return str(pins[0].with_suffix(""))
            psqs = list(p.glob("*.psq"))
            if psqs:
                return str(psqs[0].with_suffix(""))

        # Check if it's already a prefix (e.g., /path/to/db/dog_VH_FR)
        # by checking if the .pin file exists
        pin_file = Path(str(p) + ".pin")
        if pin_file.exists():
            return str(p)

        # Check if the .psq file exists
        psq_file = Path(str(p) + ".psq")
        if psq_file.exists():
            return str(p)

        # Check if it's a file (e.g., pointing to .pin directly)
        if p.is_file():
            return str(p)

        return str(p)

    def score_sequence(
        self,
        sequence: str,
        numbered: Optional[NumberedSequence] = None,
        mode: str = "framework",
    ) -> float:
        """
        Score a sequence against a target species database.

        Args:
            sequence: Amino acid sequence
            numbered: Optional NumberedSequence (will be created if None)
            mode: 'framework' for FR-only scoring, 'full' for full VH scoring

        Returns:
            T20-like score (average percent identity of top hits)
        """
        if self.blastdb_path is None:
            raise ValueError(
                "No BLAST database configured. "
                "Set blastdb_path or run 'panig download' first."
            )

        # Extract the query sequence based on mode
        if mode == "framework":
            if numbered is None:
                numbered = self._auto_number(sequence)
            if numbered is not None:
                query_seq = numbered.get_framework_sequence()
            else:
                query_seq = sequence
        else:
            query_seq = sequence

        if not query_seq:
            logger.warning("Empty query sequence")
            return 0.0

        # Run BLAST
        hits = self._run_blast(query_seq)

        # Calculate T20 score
        if not hits:
            return 0.0

        return self._calculate_t20(hits)

    def score_both(
        self,
        original: str,
        xenotypized: str,
        original_numbered: Optional[NumberedSequence] = None,
        xenotypized_numbered: Optional[NumberedSequence] = None,
        mode: str = "framework",
    ) -> Tuple[float, float]:
        """
        Score both original and animalized sequences.

        Args:
            original: Original sequence
            xenotypized: Xenotypized sequence
            original_numbered: Numbered original sequence
            xenotypized_numbered: Numbered xenotypized sequence
            mode: 'framework' or 'full'

        Returns:
            Tuple of (original_score, xenotypized_score)
        """
        # Auto-number if not provided
        if mode == "framework":
            if original_numbered is None:
                original_numbered = self._auto_number(original)
            if xenotypized_numbered is None:
                xenotypized_numbered = self._auto_number(xenotypized)

        original_score = self.score_sequence(
            original, original_numbered, mode
        )
        xenotypized_score = self.score_sequence(
            xenotypized, xenotypized_numbered, mode
        )
        return original_score, xenotypized_score

    def compute_position_coverage(
        self,
        numbered: NumberedSequence,
        species_profile,
        chain_type: str = "heavy",
    ) -> Tuple[float, int, int]:
        """
        Compute position coverage for a numbered sequence against a species profile.

        Position coverage = fraction of FR positions where the residue exists
        in the species profile (freq > 0).

        Args:
            numbered: NumberedSequence object
            species_profile: SpeciesProfile object
            chain_type: Chain type ('heavy', 'light', 'nanobody')

        Returns:
            Tuple of (coverage_fraction, covered_count, total_fr_count)
        """
        if species_profile is None:
            return 0.0, 0, 0

        total_fr = 0
        covered = 0

        for pos_info in numbered.positions:
            # Skip CDR positions
            if hasattr(pos_info, 'region') and pos_info.region.startswith("CDR"):
                continue
            # Skip if region attribute not available - use position-based CDR check
            if not hasattr(pos_info, 'region'):
                imgt_pos = pos_info.position if hasattr(pos_info, 'position') else pos_info
                if (27 <= imgt_pos <= 38) or (56 <= imgt_pos <= 65) or (105 <= imgt_pos <= 117):
                    continue

            total_fr += 1
            freq = species_profile.get_frequency(pos_info.position, pos_info.residue)
            if freq > 0:
                covered += 1

        coverage = covered / total_fr if total_fr > 0 else 0.0
        return coverage, covered, total_fr

    def _auto_number(self, sequence: str) -> Optional[NumberedSequence]:
        """Auto-number a sequence using the Numberer."""
        try:
            from panig.numbering import Numberer
            numberer = Numberer(scheme="imgt")
            return numberer.number_sequence(sequence, "temp", None)
        except Exception:
            return None

    def _run_blast(self, query_seq: str) -> List[Dict]:
        """
        Run BLASTP against the target species database.

        Args:
            query_seq: Query amino acid sequence

        Returns:
            List of hit dictionaries with 'pident', 'bitscore', etc.
        """
        # Write query to temp file
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as f:
            f.write(f">query\n{query_seq}\n")
            query_file = f.name

        # Run BLAST
        try:
            import shutil
            blastp_path = shutil.which("blastp")
            if blastp_path is None:
                # Try common conda locations
                for env_path in [
                    Path.home() / ".conda/envs/panig/bin/blastp",
                    Path.home() / ".conda/envs/llamanade/bin/blastp",
                ]:
                    if env_path.exists():
                        blastp_path = str(env_path)
                        break
            if blastp_path is None:
                raise FileNotFoundError(
                    "blastp not found on PATH. "
                    "Install with: conda install -c bioconda blast"
                )
            cmd = [
                blastp_path,
                "-query", query_file,
                "-db", self.blastdb_path,
                "-evalue", str(self.evalue),
                "-num_threads", str(self.num_threads),
                "-outfmt", "7",  # Tabular with comments
                "-max_target_seqs", "20",  # Top 20 hits
            ]

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
            )

            return self._parse_blast_output(result.stdout)

        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST failed: {e.stderr}")
            return []
        except FileNotFoundError:
            logger.error(
                "BLAST+ not found. Install with: "
                "conda install -c bioconda blast"
            )
            return []
        finally:
            # Clean up temp file
            Path(query_file).unlink(missing_ok=True)

    def _parse_blast_output(self, output: str) -> List[Dict]:
        """
        Parse BLAST tabular output (format 7).

        Args:
            output: Raw BLAST output string

        Returns:
            List of hit dictionaries
        """
        hits = []

        for line in output.strip().split("\n"):
            # Skip comment lines
            if line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) >= 12:
                try:
                    hits.append({
                        "subject_id": fields[1],
                        "pident": float(fields[2]),
                        "length": int(fields[3]),
                        "mismatch": int(fields[4]),
                        "gapopen": int(fields[5]),
                        "qstart": int(fields[6]),
                        "qend": int(fields[7]),
                        "sstart": int(fields[8]),
                        "send": int(fields[9]),
                        "evalue": float(fields[10]),
                        "bitscore": float(fields[11]),
                    })
                except (ValueError, IndexError):
                    continue

        return hits

    def _calculate_t20(self, hits: List[Dict]) -> float:
        """
        Calculate T20 score from BLAST hits.

        T20 = average percent identity of top 20 hits.

        Args:
            hits: List of BLAST hit dictionaries

        Returns:
            T20 score (0.0 to 100.0)
        """
        if not hits:
            return 0.0

        # Sort by bitscore (descending)
        sorted_hits = sorted(
            hits, key=lambda x: x["bitscore"], reverse=True
        )

        # Take top 20 hits
        top_hits = sorted_hits[:20]

        # Average percent identity
        avg_identity = sum(h["pident"] for h in top_hits) / len(top_hits)

        return avg_identity

    def download_database(
        self,
        species: str,
        chain_type: str = "VH",
        force: bool = False,
    ):
        """
        Download a species-specific BLAST database from Google Drive.

        Args:
            species: Species name
            chain_type: Chain type ('VH' or 'VHH')
            force: Force re-download even if exists
        """
        import subprocess

        cache_dir = Path.home() / ".panig" / "cache" / "blastdb"
        cache_dir.mkdir(parents=True, exist_ok=True)

        db_name = f"{species}_{chain_type}_blastdb"
        remote_path = f"gdrive:PanIg_databases/blastdb/{db_name}"

        # Check if already exists
        local_path = cache_dir / db_name
        if local_path.exists() and not force:
            logger.info(f"BLAST database already exists: {local_path}")
            self.blastdb_path = self._resolve_blastdb_path(str(local_path))
            return

        # Download
        try:
            subprocess.run(
                ["rclone", "copy", remote_path, str(cache_dir)],
                check=True,
                capture_output=True,
            )
            self.blastdb_path = str(local_path)
            logger.info(f"Downloaded BLAST database: {local_path}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to download BLAST database: {e}")
            raise

    def to_csv(
        self,
        original_score: float,
        xenotypized_score: float,
        output_path: str,
        sequence_name: str = "query",
    ):
        """
        Export scoring results to CSV.

        Args:
            original_score: Score of the original sequence
            xenotypized_score: Score of the xenotypized sequence
            output_path: Path to output CSV file
            sequence_name: Name of the sequence
        """
        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "Sequence", "Original_T20", "Xenotypized_T20", "Improvement"
            ])
            writer.writerow([
                sequence_name,
                f"{original_score:.2f}",
                f"{xenotypized_score:.2f}",
                f"{xenotypized_score - original_score:.2f}",
            ])
