"""Pytest configuration and fixtures for PanIg tests."""

import pytest
from pathlib import Path


@pytest.fixture
def test_data_dir():
    """Get the test data directory path."""
    return Path(__file__).parent / "test_data"


@pytest.fixture
def sample_nanobody_sequence():
    """Get a sample nanobody sequence for testing."""
    return (
        "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVS"
        "AISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK"
        "DFYDGSGYWGQGTQVTVSS"
    )


@pytest.fixture
def sample_antibody_sequence():
    """Get a sample antibody heavy chain sequence for testing."""
    return (
        "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVS"
        "AISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK"
        "DFYDGSGYWGQGTQVTVSS"
    )


@pytest.fixture
def sample_species_profile():
    """Get a sample species profile for testing."""
    from panig.species_profiles import SpeciesProfile

    return SpeciesProfile(
        species="dog",
        chain_type="VH",
        profile={
            1: {"E": 0.75, "D": 0.15, "Q": 0.05, "A": 0.03, "V": 0.02},
            2: {"V": 0.85, "L": 0.10, "I": 0.03, "A": 0.02},
            3: {"Q": 0.70, "E": 0.15, "K": 0.05, "R": 0.05, "D": 0.03, "N": 0.02},
        }
    )


@pytest.fixture
def tmp_fasta(tmp_path, sample_nanobody_sequence):
    """Create a temporary FASTA file with a sample sequence."""
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(f">test_nb\n{sample_nanobody_sequence}\n")
    return fasta_file


@pytest.fixture
def tmp_profile(tmp_path, sample_species_profile):
    """Create a temporary profile file."""
    profile_file = tmp_path / "test_profile.json"
    sample_species_profile.save(str(profile_file))
    return profile_file
