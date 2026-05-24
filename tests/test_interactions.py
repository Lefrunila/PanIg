"""Tests for interaction-aware substitution compatibility."""

import pytest
from panig.interactions import (
    get_compatible_residues,
    InteractionDetector,
    AA_PROPERTIES,
    IONIC_RESIDUES,
    CATION_RESIDUES,
    PI_RESIDUES,
)


class TestAAProperties:
    """Test amino acid property classification."""

    def test_all_20_aas_classified(self):
        """All 20 standard amino acids must be in AA_PROPERTIES."""
        standard = set("ACDEFGHIKLMNPQRSTVWY")
        assert set(AA_PROPERTIES.keys()) == standard

    def test_charged_positive(self):
        assert AA_PROPERTIES["K"] == "charged_positive"
        assert AA_PROPERTIES["R"] == "charged_positive"
        assert AA_PROPERTIES["H"] == "charged_positive"

    def test_charged_negative(self):
        assert AA_PROPERTIES["D"] == "charged_negative"
        assert AA_PROPERTIES["E"] == "charged_negative"

    def test_aromatic(self):
        assert AA_PROPERTIES["F"] == "aromatic"
        assert AA_PROPERTIES["W"] == "aromatic"
        assert AA_PROPERTIES["Y"] == "aromatic"


class TestGetCompatibleResidues:
    """Test the compatibility resolver."""

    # Ionic: negative residues
    def test_ionic_negative_D(self):
        result = get_compatible_residues("D", "ionic")
        assert result == {"D", "E"}

    def test_ionic_negative_E(self):
        result = get_compatible_residues("E", "ionic")
        assert result == {"D", "E"}

    # Ionic: positive residues
    def test_ionic_positive_K(self):
        result = get_compatible_residues("K", "ionic")
        assert result == {"K", "R", "H"}

    def test_ionic_positive_R(self):
        result = get_compatible_residues("R", "ionic")
        assert result == {"K", "R", "H"}

    # Cation-pi
    def test_cation_pi_cation_K(self):
        result = get_compatible_residues("K", "cation_pi")
        assert result == CATION_RESIDUES | PI_RESIDUES

    def test_cation_pi_pi_F(self):
        result = get_compatible_residues("F", "cation_pi")
        assert result == CATION_RESIDUES | PI_RESIDUES

    # Pi-pi
    def test_pi_pi_W(self):
        result = get_compatible_residues("W", "pi_pi")
        assert result == PI_RESIDUES

    def test_pi_pi_Y(self):
        result = get_compatible_residues("Y", "pi_pi")
        assert result == PI_RESIDUES

    # Unknown interaction type
    def test_unknown_type_returns_original_only(self):
        result = get_compatible_residues("A", "hydrophobic_interaction")
        assert result == {"A"}

    # Non-interacting residue in an interaction context (defensive)
    def test_non_ionic_residue_in_ionic_context(self):
        result = get_compatible_residues("A", "ionic")
        assert result == {"A"}


class TestGetInteractionMap:
    """Test the interaction map builder."""

    def test_single_type(self):
        interactions = {"ionic": {10, 25}, "cation_pi": set(), "pi_pi": set()}
        result = InteractionDetector.get_interaction_map(interactions)
        assert result == {10: "ionic", 25: "ionic"}

    def test_priority_ionic_over_cation_pi(self):
        """Position in both ionic and cation_pi should get ionic (most constrained)."""
        interactions = {"ionic": {10}, "cation_pi": {10, 30}, "pi_pi": set()}
        result = InteractionDetector.get_interaction_map(interactions)
        assert result[10] == "ionic"

    def test_priority_cation_pi_over_pi_pi(self):
        interactions = {"ionic": set(), "cation_pi": {20}, "pi_pi": {20, 40}}
        result = InteractionDetector.get_interaction_map(interactions)
        assert result[20] == "cation_pi"
        assert result[40] == "pi_pi"

    def test_empty_interactions(self):
        interactions = {"ionic": set(), "cation_pi": set(), "pi_pi": set()}
        result = InteractionDetector.get_interaction_map(interactions)
        assert result == {}

    def test_all_types(self):
        interactions = {"ionic": {1}, "cation_pi": {2, 3}, "pi_pi": {3, 4}}
        result = InteractionDetector.get_interaction_map(interactions)
        assert result == {1: "ionic", 2: "cation_pi", 3: "cation_pi", 4: "pi_pi"}
