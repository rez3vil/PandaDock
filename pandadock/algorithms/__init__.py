"""
Docking algorithms for PandaDock.

This module contains all search algorithms with a unified interface.
"""

from .base_algorithm import BaseAlgorithm
from .algorithm_factory import AlgorithmFactory

# Import specific algorithms
try:
    from .genetic_algorithm import GeneticAlgorithm, SimpleGeneticAlgorithm
    from .random_search import RandomSearchAlgorithm, SimpleRandomSearch
    from .monte_carlo import MonteCarloAlgorithm, MonteCarloSampling
    from .pandadock_algorithm import PANDADOCKAlgorithm
    from .mmff_minimization import MMFFMinimization, SimpleMMFFMinimization
    from .metal_docking import (MetalDockingAlgorithm, MetalDockingScorer, 
                               MetalCenter, MetalConstraint, MetalDockingPreparation)
    
    __all__ = [
        'BaseAlgorithm', 'AlgorithmFactory',
        'GeneticAlgorithm', 'SimpleGeneticAlgorithm',
        'RandomSearchAlgorithm', 'SimpleRandomSearch',
        'MonteCarloAlgorithm', 'MonteCarloSampling',
        'PANDADOCKAlgorithm',
        'MMFFMinimization', 'SimpleMMFFMinimization',
        'MetalDockingAlgorithm', 'MetalDockingScorer', 
        'MetalCenter', 'MetalConstraint', 'MetalDockingPreparation'
    ]
    
except ImportError as e:
    import logging
    logging.getLogger(__name__).warning(f"Some algorithms not available: {e}")
    __all__ = ['BaseAlgorithm', 'AlgorithmFactory']