"""
ANARCII wrapper for antibody/nanobody sequence numbering.

Uses ANARCII (OPIG, Oxford) for fast, accurate sequence numbering
with support for multiple schemes and chain types.

Also integrates old ANARCI for species detection and germline assignment.
"""

import logging
import subprocess
import tempfile
from typing import Dict, List, Optional

from panig.scheme import classify_position, get_scheme
from panig.sequence import NumberedPosition, NumberedSequence

logger = logging.getLogger(__name__)


class Numberer:
    """
    Antibody/nanobody sequence numberer using ANARCII + ANARCI.

    Supports:
    - Heavy chain (VH) antibodies
    - Nanobody (VHH) single-domain antibodies
    - Light chain (VL) antibodies
    - Multiple numbering schemes (IMGT, Kabat, Chothia, Martin, AHo)
    - Species detection and germline assignment via ANARCI
    """

    # ANARCII chain type mapping
    CHAIN_TYPE_MAP = {
        "H": "heavy",
        "K": "light",  # Kappa light chain
        "L": "light",  # Lambda light chain
        "A": "tcr",    # TCR alpha
        "B": "tcr",    # TCR beta
        "G": "tcr",    # TCR gamma
        "D": "tcr",    # TCR delta
    }

    def __init__(self, scheme: str = "imgt", verbose: bool = False, use_anarci: bool = True):
        """
        Initialize the numberer.

        Args:
            scheme: Numbering scheme to use ('imgt', 'kabat', 'chothia', 'martin', 'aho')
            verbose: Whether to print verbose output
            use_anarci: Whether to use old ANARCI for species detection (default: True)
        """
        self.scheme = scheme.lower()
        self.verbose = verbose
        self.scheme_def = get_scheme(self.scheme)
        self.use_anarci = use_anarci
        self._anarcii = None

    def _get_anarcii(self):
        """Lazy-load ANARCII to avoid import errors if not installed."""
        if self._anarcii is None:
            try:
                from anarcii import Anarcii
                self._anarcii = Anarcii(seq_type='antibody', mode='accuracy', cpu=True)
            except ImportError:
                raise ImportError(
                    "ANARCII is not installed. Install with: pip install anarcii"
                )
        return self._anarcii

    def _detect_species_anarci(self, sequence: str) -> dict:
        """
        Use old ANARCI to detect species and assign germlines.

        The old ANARCI (v1.3) has --assign_germline which detects species
        (human, mouse, rat, rabbit, rhesus, pig, alpaca, cow) and assigns
        V/J germline genes with identity scores.

        Args:
            sequence: Amino acid sequence

        Returns:
            Dictionary with species, v_gene, j_gene, v_identity, j_identity
        """
        result = {
            "species": None,
            "v_gene": None,
            "j_gene": None,
            "v_identity": None,
            "j_identity": None,
        }

        try:
            import anarci as old_anarci

            # Run old ANARCI with germline assignment
            # API: anarci.anarci([(name, seq)], scheme=..., output=True, assign_germline=True)
            # Returns: (numbering_output, alignment_details, hit_tables)
            # where alignment_details is [[{detail_dict}, ...], ...]
            anarci_result = old_anarci.anarci(
                [("query", sequence)],
                scheme=self.scheme,
                output=True,
                assign_germline=True,
            )

            if anarci_result and len(anarci_result) >= 2:
                alignment_details = anarci_result[1]

                if alignment_details and len(alignment_details) > 0:
                    details_list = alignment_details[0]  # First sequence's details

                    if details_list and len(details_list) > 0:
                        detail = details_list[0]  # First domain

                        # Extract species from HMM hit
                        if 'species' in detail:
                            result['species'] = detail['species']

                        # Extract germline info
                        if 'germlines' in detail:
                            germlines = detail['germlines']

                            # V gene: [('species', 'gene_name'), identity_float]
                            if 'v_gene' in germlines:
                                v_info = germlines['v_gene']
                                if isinstance(v_info, (list, tuple)) and len(v_info) >= 2:
                                    v_gene_info, v_identity = v_info
                                    if isinstance(v_gene_info, (list, tuple)) and len(v_gene_info) >= 2:
                                        result['v_gene'] = f"{v_gene_info[0]}:{v_gene_info[1]}"
                                    result['v_identity'] = v_identity

                            # J gene: [('species', 'gene_name'), identity_float]
                            if 'j_gene' in germlines:
                                j_info = germlines['j_gene']
                                if isinstance(j_info, (list, tuple)) and len(j_info) >= 2:
                                    j_gene_info, j_identity = j_info
                                    if isinstance(j_gene_info, (list, tuple)) and len(j_gene_info) >= 2:
                                        result['j_gene'] = f"{j_gene_info[0]}:{j_gene_info[1]}"
                                    result['j_identity'] = j_identity

        except ImportError:
            logger.debug("Old ANARCI not installed. Species detection disabled.")
        except Exception as e:
            logger.debug(f"ANARCI species detection failed: {e}")

        return result

    def number_sequence(
        self,
        sequence: str,
        name: str = "unknown",
        chain_type: Optional[str] = None,
    ) -> NumberedSequence:
        """
        Number a single antibody/nanobody sequence.

        Args:
            sequence: Amino acid sequence string
            name: Sequence identifier
            chain_type: Expected chain type ('heavy', 'nanobody', 'light').
                       If None, will be auto-detected by ANARCII.

        Returns:
            NumberedSequence object with numbered positions and annotations
        """
        anarcii = self._get_anarcii()

        # Run ANARCII numbering
        results = anarcii.number([sequence])

        if not results or len(results) == 0:
            raise ValueError(f"ANARCII failed to number sequence: {name}")

        # Get the first result
        result = list(results.values())[0]

        # Extract numbering and metadata
        numbering = result.get("numbering", [])
        detected_chain = result.get("chain_type", "H")
        score = result.get("score", 0.0)

        # Map chain type
        if chain_type is None:
            if detected_chain in self.CHAIN_TYPE_MAP:
                chain_type = self.CHAIN_TYPE_MAP[detected_chain]
            else:
                chain_type = "heavy"  # Default

        # For nanobodies, ANARCII may report as "H" but we need to detect VHH
        # Nanobodies typically have specific germline patterns
        if chain_type == "heavy" and self._is_nanobody(sequence, result):
            chain_type = "nanobody"

        # Build numbered positions
        positions = []
        for pos_data in numbering:
            # ANARCII v2 format: ((position, insertion), residue)
            if isinstance(pos_data, (list, tuple)) and len(pos_data) >= 2:
                pos_info, residue = pos_data[0], pos_data[1]
                if isinstance(pos_info, (list, tuple)) and len(pos_info) >= 1:
                    pos = pos_info[0]
                else:
                    pos = pos_info

                # Skip gaps
                if residue == "-" or residue == " ":
                    continue

                # Classify position into FR/CDR region
                region = classify_position(pos, self.scheme_def)

                positions.append(NumberedPosition(
                    position=pos,
                    residue=residue,
                    region=region,
                ))

        # Use ANARCI for species detection if enabled
        species = None
        v_gene = None
        j_gene = None
        v_identity = None

        if self.use_anarci:
            anarci_result = self._detect_species_anarci(sequence)
            species = anarci_result.get("species")
            v_gene = anarci_result.get("v_gene")
            j_gene = anarci_result.get("j_gene")
            v_identity = anarci_result.get("v_identity")

            if species and self.verbose:
                logger.info(f"ANARCI detected species: {species}")

        return NumberedSequence(
            name=name,
            sequence=sequence,
            chain_type=chain_type,
            scheme=self.scheme,
            positions=positions,
            species=species,
            v_gene=v_gene,
            j_gene=j_gene,
        )

    def number_sequences(
        self,
        sequences: Dict[str, str],
        chain_type: Optional[str] = None,
    ) -> List[NumberedSequence]:
        """
        Number multiple sequences in batch.

        Args:
            sequences: Dictionary of {name: sequence} pairs
            chain_type: Expected chain type for all sequences

        Returns:
            List of NumberedSequence objects
        """
        results = []
        for name, seq in sequences.items():
            try:
                numbered = self.number_sequence(seq, name, chain_type)
                results.append(numbered)
            except Exception as e:
                logger.warning(f"Failed to number {name}: {e}")
                continue
        return results

    def number_from_fasta(
        self,
        fasta_path: str,
        chain_type: Optional[str] = None,
    ) -> List[NumberedSequence]:
        """
        Number sequences from a FASTA file.

        Args:
            fasta_path: Path to FASTA file
            chain_type: Expected chain type

        Returns:
            List of NumberedSequence objects
        """
        sequences = self._parse_fasta(fasta_path)
        return self.number_sequences(sequences, chain_type)

    def _is_nanobody(self, sequence: str, anarcii_result: dict) -> bool:
        """
        Detect if a sequence is likely a nanobody (VHH).

        Nanobodies have characteristic features:
        - Specific germline genes (IGHV1 family in camelids)
        - Lack of light chain pairing residues
        - Specific CDR3 features
        """
        # Check germline assignment
        v_gene = anarcii_result.get("v_gene", "")
        if v_gene and "IGHV1" in str(v_gene):
            # Could be nanobody - check additional features
            species = anarcii_result.get("species", "")
            if species and any(s in str(species).lower() for s in [
                "alpaca", "llama", "camel", "vicuna", "guanaco"
            ]):
                return True

        # Check sequence length (nanobodies are typically 110-130 residues)
        if 100 <= len(sequence) <= 140:
            # Check for characteristic nanobody residues
            # Nanobodies have specific substitutions at key positions
            # This is a simplified heuristic
            pass

        return False

    @staticmethod
    def _parse_fasta(fasta_path: str) -> Dict[str, str]:
        """Parse a FASTA file into a dictionary of {name: sequence}."""
        sequences = {}
        current_name = None
        current_seq = []

        with open(fasta_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_name is not None:
                        sequences[current_name] = "".join(current_seq)
                    current_name = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)

        if current_name is not None:
            sequences[current_name] = "".join(current_seq)

        return sequences
