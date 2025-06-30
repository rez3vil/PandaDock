"""
Molecular structure handling for PandaDock.

This module provides unified interfaces for protein and ligand structure
operations, preparation, and manipulation.
"""

from .protein_handler import ProteinHandler, ProteinStructure
from .ligand_handler import LigandHandler
from .structure_preparation import StructurePreparation

try:
    from .flexible_residue_detector import FlexibleResidueDetector
    __all__ = [
        'ProteinHandler', 'ProteinStructure', 'LigandHandler', 'StructurePreparation',
        'FlexibleResidueDetector'
    ]
except ImportError:
    __all__ = ['ProteinHandler', 'ProteinStructure', 'LigandHandler', 'StructurePreparation']