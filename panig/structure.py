"""
Structure prediction using ImmuneBuilder.

Provides 3D structure prediction for antibodies and nanobodies
using OPIG's ImmuneBuilder (ABodyBuilder2 + NanoBodyBuilder2).
"""

import logging
import tempfile
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


class StructurePredictor:
    """
    Antibody/nanobody structure predictor using ImmuneBuilder.

    Supports:
    - ABodyBuilder2 for conventional antibodies (VH+VL)
    - NanoBodyBuilder2 for nanobodies (VHH)
    """

    def __init__(self, device: Optional[str] = None):
        """
        Initialize the structure predictor.

        Args:
            device: PyTorch device ('cpu', 'cuda', 'cuda:0', etc.)
                   If None, will auto-detect.
        """
        self.device = device
        self._abody_builder = None
        self._nanobody_builder = None

    def _get_abody_builder(self):
        """Lazy-load ABodyBuilder2."""
        if self._abody_builder is None:
            try:
                from ImmuneBuilder.ABodyBuilder2 import ABodyBuilder2
                self._abody_builder = ABodyBuilder2()
            except ImportError:
                raise ImportError(
                    "ImmuneBuilder is not installed. "
                    "Install with: pip install immunebuilder"
                )
        return self._abody_builder

    def _get_nanobody_builder(self):
        """Lazy-load NanoBodyBuilder2."""
        if self._nanobody_builder is None:
            try:
                from ImmuneBuilder.NanoBodyBuilder2 import NanoBodyBuilder2
                self._nanobody_builder = NanoBodyBuilder2()
            except ImportError:
                raise ImportError(
                    "ImmuneBuilder is not installed. "
                    "Install with: pip install immunebuilder"
                )
        return self._nanobody_builder

    def predict(
        self,
        sequence: str,
        chain_type: str = "nanobody",
        output_path: Optional[str] = None,
        name: str = "predicted",
    ) -> str:
        """
        Predict 3D structure for an antibody or nanobody.

        Args:
            sequence: Amino acid sequence
            chain_type: 'nanobody' for VHH, 'antibody' for VH+VL
            output_path: Path to save PDB file. If None, uses temp directory.
            name: Name for the output file

        Returns:
            Path to the predicted PDB file
        """
        if chain_type == "nanobody":
            return self._predict_nanobody(sequence, output_path, name)
        elif chain_type == "antibody":
            return self._predict_antibody(sequence, output_path, name)
        else:
            raise ValueError(
                f"Unsupported chain_type: {chain_type}. "
                f"Use 'nanobody' or 'antibody'."
            )

    def _predict_nanobody(
        self,
        sequence: str,
        output_path: Optional[str],
        name: str,
    ) -> str:
        """Predict nanobody structure using NanoBodyBuilder2."""
        builder = self._get_nanobody_builder()

        if output_path is None:
            output_path = tempfile.mkdtemp()

        Path(output_path).mkdir(parents=True, exist_ok=True)
        output_file = Path(output_path) / f"{name}.pdb"

        # NanoBodyBuilder2 expects a dictionary with 'H' key for heavy chain
        sequences = {"H": sequence}
        nanobody = builder.predict(sequences)
        nanobody.save(str(output_file))

        logger.info(f"Predicted nanobody structure: {output_file}")
        return str(output_file)

    def _predict_antibody(
        self,
        sequence: str,
        output_path: Optional[str],
        name: str,
    ) -> str:
        """
        Predict antibody structure using ABodyBuilder2.

        Note: For conventional antibodies, both heavy and light chain
        sequences are needed. If only heavy chain is provided, it will
        predict the heavy chain Fv region only.
        """
        builder = self._get_abody_builder()

        if output_path is None:
            output_path = tempfile.mkdtemp()

        Path(output_path).mkdir(parents=True, exist_ok=True)
        output_file = Path(output_path) / f"{name}.pdb"

        # ABodyBuilder2 expects a dictionary with 'H' and optionally 'L' keys
        sequences = {"H": sequence}
        antibody = builder.predict(sequences)
        antibody.save(str(output_file))

        logger.info(f"Predicted antibody structure: {output_file}")
        return str(output_file)

    def predict_with_confidence(
        self,
        sequence: str,
        chain_type: str = "nanobody",
        output_path: Optional[str] = None,
        name: str = "predicted",
    ) -> dict:
        """
        Predict structure with confidence estimates.

        Args:
            sequence: Amino acid sequence
            chain_type: 'nanobody' or 'antibody'
            output_path: Path to save PDB file
            name: Name for the output file

        Returns:
            Dictionary with 'pdb_path' and 'confidence' scores
        """
        pdb_path = self.predict(sequence, chain_type, output_path, name)

        # Parse confidence from PDB B-factors
        confidence = self._extract_confidence(pdb_path)

        return {
            "pdb_path": pdb_path,
            "confidence": confidence,
        }

    def _extract_confidence(self, pdb_path: str) -> dict:
        """
        Extract confidence scores from PDB B-factors.

        ImmuneBuilder stores confidence estimates as B-factors.
        Lower B-factor = higher confidence.
        """
        residue_confidence = {}

        try:
            with open(pdb_path, "r") as f:
                for line in f:
                    if line.startswith("ATOM") and line[12:16].strip() == "CA":
                        res_num = int(line[22:26].strip())
                        b_factor = float(line[60:66].strip())
                        residue_confidence[res_num] = b_factor
        except Exception as e:
            logger.warning(f"Could not parse confidence from PDB: {e}")

        return residue_confidence

    def predict_ensemble(
        self,
        sequence: str,
        chain_type: str = "nanobody",
        n_models: int = 5,
        output_path: Optional[str] = None,
        name: str = "ensemble",
    ) -> list:
        """
        Predict an ensemble of structures for uncertainty estimation.

        Args:
            sequence: Amino acid sequence
            chain_type: 'nanobody' or 'antibody'
            n_models: Number of models in the ensemble
            output_path: Path to save PDB files
            name: Base name for output files

        Returns:
            List of PDB file paths
        """
        if output_path is None:
            output_path = tempfile.mkdtemp()

        pdb_files = []
        for i in range(n_models):
            model_name = f"{name}_model_{i+1}"
            pdb_path = self.predict(
                sequence, chain_type, output_path, model_name
            )
            pdb_files.append(pdb_path)

        logger.info(
            f"Generated {n_models} structure models in {output_path}"
        )
        return pdb_files
