"""
PandaDock: Modular, Multi-Strategy, High-Performance Molecular Docking Software

A comprehensive molecular docking framework supporting multiple docking strategies:
- Precise: Physics-based docking (Glide-style)
- Balanced: ML-based docking (DiffDock-style) 
- Fast: GA-based docking (Vina-style)

Features:
- Flexible protein-ligand docking with side-chain flexibility
- Multiple scoring functions and ML rescoring
- Comprehensive pose filtering and validation
- Multi-format support (SMILES, SDF, MOL2, PDB)
- HTML reporting with pose visualization
- GPU acceleration and parallel processing
"""

__version__ = "1.0.0"
__author__ = "PandaDock Development Team"
__email__ = "info@pandadock.org"

from .config import PandaDockConfig
from .docking.base_engine import DockingEngine, Pose
from .docking.physics_engine import PhysicsEngine
from .docking.ml_engine import MLEngine  
from .docking.ga_engine import GAEngine

__all__ = [
    'PandaDockConfig',
    'DockingEngine',
    'Pose',
    'PhysicsEngine',
    'MLEngine', 
    'GAEngine',
    '__version__',
]