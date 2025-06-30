"""
PandaDock - Refactored Architecture

This is the new modular, maintainable version of PandaDock with:
- Clean separation of concerns
- CPU/GPU transparency 
- Robust error handling
- Easy debugging and extensibility
"""

__version__ = "2.0.0-refactored"

# Core modules
from .core import DockingEngine, PoseGenerator, ResultProcessor

# Hardware abstraction
from .hardware import DeviceManager, ComputeBackend, PerformanceMonitor

# Algorithm system
from .algorithms import BaseAlgorithm, AlgorithmFactory

# Scoring system  
from .scoring import BaseScoringFunction, ScoringFunctionFactory

# Molecule handling
from .molecules import ProteinHandler, LigandHandler, StructurePreparation

# I/O operations
from .io import ResultWriters, ReportGenerators

# CLI interface
from .cli import create_argument_parser, get_config_from_args

# Virtual screening and batch processing
from .screening import VirtualScreeningManager, BatchScreeningManager

__all__ = [
    # Version
    '__version__',
    
    # Core components
    'DockingEngine', 'PoseGenerator', 'ResultProcessor',
    
    # Hardware abstraction
    'DeviceManager', 'ComputeBackend', 'PerformanceMonitor',
    
    # Algorithm system
    'BaseAlgorithm', 'AlgorithmFactory',
    
    # Scoring system
    'BaseScoringFunction', 'ScoringFunctionFactory',
    
    # Molecule handling
    'ProteinHandler', 'LigandHandler', 'StructurePreparation',
    
    # I/O operations
    'ResultWriters', 'ReportGenerators',
    
    # CLI interface
    'create_argument_parser', 'get_config_from_args',
    
    # Virtual screening and batch processing
    'VirtualScreeningManager', 'BatchScreeningManager'
]