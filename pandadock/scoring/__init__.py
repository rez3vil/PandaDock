"""
Scoring functions for PandaDock.

This module contains all scoring functions with a unified interface.
"""

from .base_scoring import BaseScoringFunction
from .scoring_factory import ScoringFunctionFactory

# Import available scoring functions
try:
    from .physics_based_scoring import PhysicsBasedScoringFunction
    PHYSICS_AVAILABLE = True
except ImportError:
    PHYSICS_AVAILABLE = False

# Try to import unified scoring functions
try:
    from .unified_scoring import CompositeScoringFunction, EnhancedScoringFunction
    UNIFIED_AVAILABLE = True
    COMPOSITE_AVAILABLE = True
    ENHANCED_AVAILABLE = True
except ImportError:
    UNIFIED_AVAILABLE = False
    COMPOSITE_AVAILABLE = False
    ENHANCED_AVAILABLE = False

try:
    from .physics_scoring import PhysicsScoringFunction
    LEGACY_PHYSICS_AVAILABLE = True
except ImportError:
    LEGACY_PHYSICS_AVAILABLE = False

# Build __all__ list based on what's available
__all__ = ['BaseScoringFunction', 'ScoringFunctionFactory']

if PHYSICS_AVAILABLE:
    __all__.append('PhysicsBasedScoringFunction')
if UNIFIED_AVAILABLE:
    __all__.extend(['CompositeScoringFunction', 'EnhancedScoringFunction'])
if LEGACY_PHYSICS_AVAILABLE:
    __all__.append('PhysicsScoringFunction')