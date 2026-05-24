"""
Structural interaction detection using Protinter.

Identifies intra-chain interactions (ionic, cation-pi, pi-pi)
that should be preserved during xenotypization.
"""

import logging
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Set, Union

logger = logging.getLogger(__name__)


# Residues involved in specific interaction types
IONIC_RESIDUES = {"D", "E", "K", "R"}  # Asp, Glu, Lys, Arg
CATION_RESIDUES = {"K", "R", "H"}  # Lys, Arg, His
PI_RESIDUES = {"F", "W", "Y", "H"}  # Phe, Trp, Tyr, His

# Chemical property classification for all 20 standard amino acids
AA_PROPERTIES = {
    "K": "charged_positive", "R": "charged_positive", "H": "charged_positive",
    "D": "charged_negative", "E": "charged_negative",
    "S": "polar", "T": "polar", "N": "polar", "Q": "polar", "C": "polar",
    "A": "hydrophobic", "V": "hydrophobic", "I": "hydrophobic",
    "L": "hydrophobic", "M": "hydrophobic", "P": "hydrophobic",
    "F": "aromatic", "W": "aromatic", "Y": "aromatic",
    "G": "special",
}


def get_compatible_residues(
    original_aa: str,
    interaction_type: str,
) -> Set[str]:
    """
    Get residues chemically compatible with original_aa for a given interaction type.

    For ionic interactions: only residues with the same charge sign.
    For cation-pi: union of cation and pi groups (can't determine partner role
    without structural context).
    For pi-pi: only other aromatic residues.

    Args:
        original_aa: The current amino acid at the interacting position.
        interaction_type: One of 'ionic', 'cation_pi', 'pi_pi'.

    Returns:
        Set of single-letter amino acid codes that are compatible replacements.
        Includes original_aa itself.
    """
    if interaction_type == "ionic":
        if original_aa in {"D", "E"}:
            return {"D", "E"}
        elif original_aa in {"K", "R", "H"}:
            return {"K", "R", "H"}
        else:
            return {original_aa}
    elif interaction_type == "cation_pi":
        return CATION_RESIDUES | PI_RESIDUES
    elif interaction_type == "pi_pi":
        return PI_RESIDUES
    else:
        return {original_aa}


class InteractionDetector:
    """
    Detect intra-chain interactions in antibody structures.

    Uses Protinter to identify:
    - Ionic interactions (salt bridges)
    - Cation-pi interactions
    - Pi-pi stacking interactions

    Residues involved in these interactions should be excluded
    from xenotypization to preserve structural stability.
    """

    def __init__(self, distance_cutoff: float = 5.0):
        """
        Initialize the interaction detector.

        Args:
            distance_cutoff: Maximum distance (Angstroms) for interaction detection
        """
        self.distance_cutoff = distance_cutoff

    def detect_interactions(
        self,
        pdb_path: Optional[str] = None,
        sequence: Optional[str] = None,
        numbered_positions: Optional[list] = None,
    ) -> Dict[str, Set[int]]:
        """
        Detect structural interactions in an antibody structure.

        Args:
            pdb_path: Path to PDB file
            sequence: Amino acid sequence (used if no PDB)
            numbered_positions: List of NumberedPosition objects

        Returns:
            Dictionary mapping interaction type to set of positions involved
            {'ionic': {10, 25}, 'cation_pi': {30}, 'pi_pi': {45, 67}}
        """
        if pdb_path:
            return self._detect_from_pdb(pdb_path)
        else:
            # Without structure, use sequence-based heuristics
            return self.detect_from_sequence(sequence, numbered_positions)

    def _detect_from_pdb(self, pdb_path: str) -> Dict[str, Set[int]]:
        """
        Detect interactions from a PDB file using Protinter.

        Protinter has no clean Python API that returns data -- its main()
        and interlib.calc_inter() only print to stdout. We use BioPython's
        PDBParser + protinter.interlib functions directly to capture the
        interacting residue pairs instead of just printing them.

        Args:
            pdb_path: Path to PDB file

        Returns:
            Dictionary of interaction types to positions
        """
        interactions = {
            "ionic": set(),
            "cation_pi": set(),
            "pi_pi": set(),
        }

        try:
            from Bio.PDB.PDBParser import PDBParser
            from protinter.interlib import get_res, calc_inter, amino, center_mass
        except ImportError:
            raise RuntimeError(
                "Protinter/BioPython not installed. "
                "Install with: pip install protinter biopython"
            )

        try:
            p = PDBParser(QUIET=True)
            structure = p.get_structure("X", pdb_path)

            interaction_types = {
                "ionic": {"amino_type": "ionic", "distmax": 6.0},
                "cation_pi": {"amino_type": "cationpi", "distmax": 6.0},
                "pi_pi": {"amino_type": "aroaro", "distmin": 4.5, "distmax": 7.0},
            }

            for model in structure:
                for chain in model:
                    for int_type, params in interaction_types.items():
                        amino_type = params["amino_type"]

                        # Compute center of mass for aromatic residues if needed
                        if amino_type in ("aroaro", "arosul"):
                            for resid in chain:
                                if resid.get_resname() in amino["aroaro"]:
                                    center_mass(resid)

                        residue_dict = get_res(chain, amino_type=amino_type)
                        if not residue_dict:
                            continue

                        found = self._calc_inter_capture(
                            residue_dict=residue_dict,
                            amino_type=amino_type,
                            distmin=params.get("distmin", 0),
                            distmax=params["distmax"],
                            interval=0,
                            atommindist=30,
                        )

                        for resid1, resid2 in found:
                            pos1 = self._residue_seq_number(resid1)
                            pos2 = self._residue_seq_number(resid2)
                            if pos1 is not None:
                                interactions[int_type].add(pos1)
                            if pos2 is not None:
                                interactions[int_type].add(pos2)

            logger.info(
                f"Protinter detected: {len(interactions['ionic'])} ionic, "
                f"{len(interactions['cation_pi'])} cation-pi, "
                f"{len(interactions['pi_pi'])} pi-pi residues"
            )

        except Exception as e:
            raise RuntimeError(f"Protinter analysis failed: {e}") from e

        return interactions

    @staticmethod
    def _residue_seq_number(resid) -> Optional[int]:
        """Extract integer sequence number from a BioPython Residue object."""
        try:
            # BioPython Residue id is a tuple like (' ', 42, ' ')
            resid_id = resid.get_id()
            if isinstance(resid_id, tuple) and len(resid_id) >= 2:
                return int(resid_id[1])
            # Fallback: parse from string repr  "Residue ARG id=42"
            return int(str(resid).split("=")[1].strip())
        except (ValueError, IndexError, AttributeError):
            return None

    @staticmethod
    def _calc_inter_capture(
        residue_dict: dict,
        amino_type: str,
        distmin: float = 0,
        distmax: float = 7.0,
        distON: float = 3.5,
        distS: float = 4.0,
        interval: int = 0,
        atommindist: float = 30,
    ) -> list:
        """
        Reproduce protinter.interlib.calc_inter logic but return the list
        of (resid1, resid2) BioPython residue pairs instead of printing.

        This is adapted from protinter v0.9.2 interlib.calc_inter.
        """
        from protinter.interlib import (
            hydrophobicfun,
            disulphidefun,
            ionicfun,
            cationpifun,
            aroarofun,
            arosulfun,
            hbond_main_mainfun,
            hbond_main_sidefun,
            hbond_side_sidefun,
            within_radiusfun,
        )

        found = []
        keys = sorted([int(x) for x in list(residue_dict.keys())])

        for i in keys:
            for j in keys:
                if i != j and abs(i - j) > interval:
                    InterOfResisFound = False
                    resid1 = residue_dict[str(i)]
                    resid2 = residue_dict[str(j)]
                    ResiduesTooFarApart = False
                    WithinMinimumDist = False
                    res = None

                    for atom1 in resid1:
                        if amino_type in (
                            "hbond_main_side",
                            "hbond_side_side",
                            "hbond_main_main",
                        ):
                            if atom1.get_name() == "C":
                                continue
                        if InterOfResisFound:
                            break
                        if ResiduesTooFarApart:
                            ResiduesTooFarApart = False
                            res = None
                            break
                        for atom2 in resid2:
                            if not WithinMinimumDist:
                                if abs(atom1 - atom2) > atommindist:
                                    ResiduesTooFarApart = True
                                    res = None
                                    break
                                else:
                                    WithinMinimumDist = True
                            if amino_type == "hydrophobic":
                                res = hydrophobicfun(atom1, atom2, dist=distmax)
                            elif amino_type == "disulphide":
                                res = disulphidefun(atom1, atom2, dist=distmax)
                            elif amino_type == "ionic":
                                res = ionicfun(atom1, atom2, resid1, resid2, dist=distmax)
                            elif amino_type == "cationpi":
                                res = cationpifun(atom1, atom2, resid1, resid2, dist=distmax)
                            elif amino_type == "hbond_main_main":
                                res = hbond_main_mainfun(atom1, atom2, distON=distON, distS=distS)
                            elif amino_type == "hbond_main_side":
                                res = hbond_main_sidefun(atom1, atom2, distON=distON, distS=distS)
                            elif amino_type == "hbond_side_side":
                                res = hbond_side_sidefun(atom1, atom2, distON=distON, distS=distS)
                            elif amino_type == "within_radius":
                                res = within_radiusfun(atom1, atom2, dist_min=distmax)
                            if res:
                                InterOfResisFound = True
                                break

                    if amino_type == "aroaro":
                        res = aroarofun(resid1, resid2, dmin=distmin, dmax=distmax)
                    elif amino_type == "arosul":
                        res = arosulfun(resid1, resid2, dist=distmax)

                    if res:
                        if (resid1, resid2) not in found and (resid2, resid1) not in found:
                            found.append((resid1, resid2))

        return found

    def detect_from_sequence(
        self,
        sequence: Optional[str],
        numbered_positions: Optional[list],
    ) -> Dict[str, Set[int]]:
        """
        Detect potential interactions from sequence alone.

        Uses heuristics based on residue types and proximity in sequence.
        Less accurate than structure-based detection but works without PDB.

        Args:
            sequence: Amino acid sequence
            numbered_positions: List of NumberedPosition objects

        Returns:
            Dictionary of interaction types to positions
        """
        interactions = {
            "ionic": set(),
            "cation_pi": set(),
            "pi_pi": set(),
        }

        if numbered_positions is None:
            return interactions

        # Group residues by type
        ionic_positions = []
        cation_positions = []
        pi_positions = []

        for pos in numbered_positions:
            if pos.residue in IONIC_RESIDUES:
                ionic_positions.append(pos.position)
            if pos.residue in CATION_RESIDUES:
                cation_positions.append(pos.position)
            if pos.residue in PI_RESIDUES:
                pi_positions.append(pos.position)

        # Heuristic: residues of complementary types within 5 positions
        # in sequence are likely to form interactions
        for i, pos1 in enumerate(ionic_positions):
            for j, pos2 in enumerate(ionic_positions):
                if i < j and abs(pos1 - pos2) <= 5:
                    # Check for complementary charges
                    aa1 = self._get_residue_at(pos1, numbered_positions)
                    aa2 = self._get_residue_at(pos2, numbered_positions)
                    if self._are_complementary_charges(aa1, aa2):
                        interactions["ionic"].add(pos1)
                        interactions["ionic"].add(pos2)

        # Cation-pi heuristic
        for cat_pos in cation_positions:
            for pi_pos in pi_positions:
                if abs(cat_pos - pi_pos) <= 5:
                    interactions["cation_pi"].add(cat_pos)
                    interactions["cation_pi"].add(pi_pos)

        # Pi-pi stacking heuristic
        for i, pos1 in enumerate(pi_positions):
            for j, pos2 in enumerate(pi_positions):
                if i < j and abs(pos1 - pos2) <= 5:
                    interactions["pi_pi"].add(pos1)
                    interactions["pi_pi"].add(pos2)

        return interactions

    def _get_residue_at(
        self, position: int, numbered_positions: list
    ) -> Optional[str]:
        """Get the amino acid at a specific position."""
        for pos in numbered_positions:
            if pos.position == position:
                return pos.residue
        return None

    def _are_complementary_charges(self, aa1: str, aa2: str) -> bool:
        """Check if two amino acids have complementary charges."""
        positive = {"K", "R", "H"}
        negative = {"D", "E"}
        return (
            (aa1 in positive and aa2 in negative) or
            (aa1 in negative and aa2 in positive)
        )

    def get_excluded_positions(
        self,
        interactions: Dict[str, Set[int]],
    ) -> Set[int]:
        """
        Get all positions that should be excluded from xenotypization.

        Args:
            interactions: Dictionary of interaction types to positions

        Returns:
            Set of all positions involved in any interaction
        """
        excluded = set()
        for positions in interactions.values():
            excluded.update(positions)
        return excluded

    @staticmethod
    def get_interaction_map(
        interactions: Dict[str, Set[int]],
    ) -> Dict[int, str]:
        """
        Get a mapping of positions to their interaction types.

        When a position is involved in multiple interaction types,
        the most constrained type is used (ionic > cation_pi > pi_pi).

        Args:
            interactions: Dictionary of interaction types to positions

        Returns:
            Dict mapping position -> interaction_type
        """
        priority = ["ionic", "cation_pi", "pi_pi"]
        result: Dict[int, str] = {}
        for int_type in priority:
            for pos in interactions.get(int_type, set()):
                if pos not in result:
                    result[pos] = int_type
        return result
