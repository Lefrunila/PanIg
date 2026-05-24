"""
PanIg: Pan-species Immunoglobulin Xenotypization Tool

A tool for adapting antibody and nanobody framework regions to match
the sequence preferences of target veterinary species.

Supports both:
- Xenotypization: Adapt sequences to a target species (inverse of humanization)
- Humanization: Adapt sequences to human preferences

Species-specific terms: caninize (dog), felinize (cat), equinize (horse), bovinize (cattle)

Licensed under MIT License.
"""

__version__ = "0.1.0"
__author__ = "PanIg Contributors"

from panig.numbering import Numberer
from panig.sequence import NumberedSequence
from panig.xenotypizer import Xenotypizer, XenotypizationResult
from panig.species_profiles import SpeciesProfile
from panig.structure import StructurePredictor
from panig.scorer import Scorer

# Backward compatibility aliases
Animalizer = Xenotypizer
AnimalizationResult = XenotypizationResult

__all__ = [
    "Numberer",
    "NumberedSequence",
    "Xenotypizer",
    "XenotypizationResult",
    # Backward compatibility
    "Animalizer",
    "AnimalizationResult",
    "SpeciesProfile",
    "StructurePredictor",
    "Scorer",
]
