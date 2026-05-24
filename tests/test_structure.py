"""Tests for the structure prediction module."""

import os
import shutil
import tempfile

import pytest

from panig.structure import StructurePredictor


# Shared VHH sequence for structure tests
VHH_SEQUENCE = (
    "QVQLVESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAITWSGGNTYY"
    "ADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAADRGYYGSGYWGQGTQVTVSS"
)

# Check if blastp is available for integration tests
_BLASTP_AVAILABLE = shutil.which("blastp") is not None


def _can_load_nanobody_builder():
    """Check if ImmuneBuilder NanoBodyBuilder2 can actually be loaded."""
    try:
        from ImmuneBuilder.NanoBodyBuilder2 import NanoBodyBuilder2
        NanoBodyBuilder2()
        return True
    except Exception:
        return False


_CAN_LOAD_BUILDER = _can_load_nanobody_builder()


class TestStructurePredictor:
    """Test StructurePredictor class."""

    def test_instantiation_default(self):
        """Test default instantiation."""
        predictor = StructurePredictor()
        assert predictor.device is None
        assert predictor._abody_builder is None
        assert predictor._nanobody_builder is None

    def test_instantiation_cpu(self):
        """Test instantiation with CPU device."""
        predictor = StructurePredictor(device="cpu")
        assert predictor.device == "cpu"

    def test_predict_invalid_chain_type(self):
        """Test predict with invalid chain type raises ValueError."""
        predictor = StructurePredictor(device="cpu")
        with pytest.raises(ValueError, match="Unsupported chain_type"):
            predictor.predict("ACDEFG", chain_type="invalid")

    def test_extract_confidence_on_sample_pdb(self, tmp_path):
        """Test _extract_confidence on a synthetic PDB file."""
        predictor = StructurePredictor()

        # Create a minimal PDB with CA atoms and B-factors
        pdb_content = (
            "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 95.50           C\n"
            "ATOM      2  CA  GLY A   2       4.000   5.000   6.000  1.00 88.30           C\n"
            "ATOM      3  CA  SER A   3       7.000   8.000   9.000  1.00 72.10           C\n"
            "END\n"
        )
        pdb_file = tmp_path / "test.pdb"
        pdb_file.write_text(pdb_content)

        confidence = predictor._extract_confidence(str(pdb_file))

        assert len(confidence) == 3
        assert confidence[1] == pytest.approx(95.50)
        assert confidence[2] == pytest.approx(88.30)
        assert confidence[3] == pytest.approx(72.10)

    def test_extract_confidence_empty_pdb(self, tmp_path):
        """Test _extract_confidence on PDB with no ATOM records."""
        predictor = StructurePredictor()

        pdb_file = tmp_path / "empty.pdb"
        pdb_file.write_text("REMARK   1 EMPTY PDB\nEND\n")

        confidence = predictor._extract_confidence(str(pdb_file))
        assert confidence == {}

    def test_extract_confidence_nonexistent_file(self, tmp_path):
        """Test _extract_confidence handles missing file gracefully."""
        predictor = StructurePredictor()

        confidence = predictor._extract_confidence(
            str(tmp_path / "nonexistent.pdb")
        )
        assert confidence == {}


@pytest.mark.skipif(
    not _CAN_LOAD_BUILDER,
    reason="ImmuneBuilder NanoBodyBuilder2 cannot be loaded (model corrupted or missing)",
)
class TestStructurePredictIntegration:
    """Integration tests for structure prediction with ImmuneBuilder."""

    @pytest.fixture(autouse=True)
    def _setup_tmpdir(self):
        self.tmpdir = tempfile.mkdtemp()
        yield

    def test_predict_nanobody_creates_pdb(self):
        """Test that predict() creates a valid PDB file for a VHH."""
        predictor = StructurePredictor(device="cpu")

        pdb_path = predictor.predict(
            VHH_SEQUENCE,
            chain_type="nanobody",
            output_path=self.tmpdir,
            name="test_vhh",
        )

        # File must exist
        assert os.path.isfile(pdb_path), f"PDB file not created: {pdb_path}"

        # Must contain ATOM records
        with open(pdb_path) as f:
            lines = f.readlines()

        atom_lines = [l for l in lines if l.startswith("ATOM")]
        assert len(atom_lines) > 0, "PDB file has no ATOM records"

        # Must contain at least one CA atom
        ca_lines = [l for l in atom_lines if l[12:16].strip() == "CA"]
        assert len(ca_lines) > 0, "PDB file has no CA atoms"

    def test_predict_with_confidence(self):
        """Test predict_with_confidence returns dict with pdb_path and confidence."""
        predictor = StructurePredictor(device="cpu")

        result = predictor.predict_with_confidence(
            VHH_SEQUENCE,
            chain_type="nanobody",
            output_path=self.tmpdir,
            name="test_vhh_conf",
        )

        assert "pdb_path" in result
        assert "confidence" in result
        assert os.path.isfile(result["pdb_path"])

        # Confidence should be a non-empty dict of residue_num -> b_factor
        conf = result["confidence"]
        assert isinstance(conf, dict)
        assert len(conf) > 0

        # All confidence values should be non-negative floats
        for res_num, b_factor in conf.items():
            assert isinstance(res_num, int)
            assert isinstance(b_factor, float)
            assert b_factor >= 0.0
