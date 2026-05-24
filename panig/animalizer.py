"""
Backward compatibility wrapper for panig.animalizer.

The canonical module is now panig.xenotypizer. This module re-exports
everything so that existing imports continue to work.
"""

from panig.xenotypizer import *  # noqa: F401,F403
from panig.xenotypizer import (  # noqa: F401
    SubstitutionRecord,
    XenotypizationResult,
    Xenotypizer,
    # Backward compat aliases
    AnimalizationResult,
    Animalizer,
)
