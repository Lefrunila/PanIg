"""
Backward compatibility wrapper for panig.vhh_animalizer.

The canonical module is now panig.vhh_xenotypizer. This module re-exports
everything so that existing imports continue to work.
"""

from panig.vhh_xenotypizer import *  # noqa: F401,F403
from panig.vhh_xenotypizer import (  # noqa: F401
    VHHSubstitutionRecord,
    VHHXenotypizationResult,
    VHHXenotypizer,
    VHH_LOCKED_POSITIONS,
    # Backward compat aliases
    VHHAnimalizationResult,
    VHHAnimalizer,
)
