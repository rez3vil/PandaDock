"""
PandaDock: Modular, Multi-Strategy, High-Performance Docking Software

A comprehensive molecular docking framework supporting:
- Physics-based docking (Glide-style)
- Deep learning-based docking (DiffDock/Boltz-style)
- Fast virtual screening (AutoDock Vina-style GA search)
- Flexible protein-ligand docking
- Clash-free pose refinement
- Comprehensive reporting with IC50, LE, ΔG
"""

__version__ = "1.0.0"
__author__ = "PandaDock Development Team"

from .docking.physics_engine import PhysicsEngine
from .docking.ml_engine import MLEngine
from .docking.ga_engine import GAEngine
from .docking.flexible_docking import FlexibleDocking
from .docking.pose_filtering import PoseFiltering
from .scoring.scoring_functions import ScoringFunctions
from .utils.ic50_calculator import IC50Calculator

__all__ = [
    'PhysicsEngine',
    'MLEngine', 
    'GAEngine',
    'FlexibleDocking',
    'PoseFiltering',
    'ScoringFunctions',
    'IC50Calculator'
]