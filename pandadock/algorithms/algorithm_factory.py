"""
Factory for creating docking algorithms.

This module provides a centralized way to create different docking
algorithms with consistent interfaces.
"""

from typing import Any, Dict, List, Optional
import logging

from ..hardware import ComputeBackend
from .base_algorithm import BaseAlgorithm


class AlgorithmFactory:
    """Factory for creating docking algorithms."""
    
    def __init__(self, compute_backend: ComputeBackend = None):
        """Initialize the algorithm factory."""
        self.compute_backend = compute_backend
        self.logger = logging.getLogger(__name__)
    
    def create_algorithm(self, algorithm_type: str, scoring_function: Any, 
                        **kwargs) -> BaseAlgorithm:
        """
        Create a docking algorithm.
        
        Args:
            algorithm_type: Type of algorithm to create
            scoring_function: Scoring function for the algorithm
            **kwargs: Algorithm-specific parameters
            
        Returns:
            Configured algorithm instance
        """
        self.logger.info(f"Creating {algorithm_type} algorithm")
        
        try:
            if algorithm_type == 'genetic':
                return self._create_genetic_algorithm(scoring_function, **kwargs)
            elif algorithm_type == 'random':
                return self._create_random_search(scoring_function, **kwargs)
            elif algorithm_type == 'monte-carlo':
                return self._create_monte_carlo(scoring_function, **kwargs)
            elif algorithm_type == 'pandadock':
                return self._create_pandadock(scoring_function, **kwargs)
            elif algorithm_type == 'metal':
                return self._create_metal_docking(scoring_function, **kwargs)
            else:
                raise ValueError(f"Unknown algorithm type: {algorithm_type}")
                
        except Exception as e:
            self.logger.error(f"Failed to create {algorithm_type} algorithm: {e}")
            # Fallback to genetic algorithm
            self.logger.info("Falling back to genetic algorithm")
            return self._create_genetic_algorithm(scoring_function, **kwargs)
    
    def _create_genetic_algorithm(self, scoring_function: Any, **kwargs) -> BaseAlgorithm:
        """Create genetic algorithm."""
        try:
            from .genetic_algorithm import GeneticAlgorithm
            return GeneticAlgorithm(scoring_function, **kwargs)
        except ImportError:
            from .genetic_algorithm import SimpleGeneticAlgorithm
            return SimpleGeneticAlgorithm(scoring_function, **kwargs)
    
    def _create_random_search(self, scoring_function: Any, **kwargs) -> BaseAlgorithm:
        """Create random search algorithm."""
        try:
            from .random_search import RandomSearchAlgorithm
            return RandomSearchAlgorithm(scoring_function, **kwargs)
        except ImportError:
            from .random_search import SimpleRandomSearch
            return SimpleRandomSearch(scoring_function, **kwargs)
    
    def _create_monte_carlo(self, scoring_function: Any, **kwargs) -> BaseAlgorithm:
        """Create Monte Carlo algorithm."""
        try:
            from .monte_carlo import MonteCarloAlgorithm
            return MonteCarloAlgorithm(scoring_function, **kwargs)
        except ImportError:
            # Fallback to genetic algorithm
            self.logger.warning("Monte Carlo not available, using genetic algorithm")
            return self._create_genetic_algorithm(scoring_function, **kwargs)
    
    def _create_pandadock(self, scoring_function: Any, **kwargs) -> BaseAlgorithm:
        """Create PANDADOCK algorithm."""
        try:
            from .pandadock_algorithm import PANDADOCKAlgorithm
            return PANDADOCKAlgorithm(scoring_function, **kwargs)
        except ImportError:
            # Fallback to genetic algorithm
            self.logger.warning("PANDADOCK not available, using genetic algorithm")
            return self._create_genetic_algorithm(scoring_function, **kwargs)
    
    def _create_metal_docking(self, scoring_function: Any, **kwargs) -> BaseAlgorithm:
        """Create metal-aware docking algorithm."""
        try:
            from .metal_docking import MetalDockingAlgorithm, MetalDockingPreparation
            
            # Extract metal-specific parameters
            protein = kwargs.get('protein')
            metal_centers = kwargs.get('metal_centers')
            metal_constraints = kwargs.get('metal_constraints')
            
            # Auto-detect metal centers if not provided
            if protein and not metal_centers:
                metal_centers = MetalDockingPreparation.detect_metal_centers(protein)
                self.logger.info(f"Auto-detected {len(metal_centers)} metal centers")
            
            # Create metal constraints if not provided
            if metal_centers and not metal_constraints:
                metal_constraints = MetalDockingPreparation.create_metal_constraints(metal_centers)
                self.logger.info(f"Created {len(metal_constraints)} metal constraints")
            
            # Create base algorithm
            base_algorithm_type = kwargs.get('base_algorithm', 'genetic')
            base_kwargs = {k: v for k, v in kwargs.items() 
                          if k not in ['protein', 'metal_centers', 'metal_constraints', 'base_algorithm']}
            base_algorithm = self.create_algorithm(base_algorithm_type, scoring_function, **base_kwargs)
            
            # Create metal docking algorithm
            return MetalDockingAlgorithm(
                base_algorithm=base_algorithm,
                metal_centers=metal_centers or [],
                metal_constraints=metal_constraints or [],
                refinement_iterations=kwargs.get('refinement_iterations', 20)
            )
            
        except ImportError:
            self.logger.warning("Metal docking not available, using genetic algorithm")
            return self._create_genetic_algorithm(scoring_function, **kwargs)
    
    def create_mmff_minimizer(self, **kwargs):
        """Create MMFF94 minimization instance."""
        try:
            from .mmff_minimization import MMFFMinimization
            return MMFFMinimization(**kwargs)
        except ImportError:
            from .mmff_minimization import SimpleMMFFMinimization
            return SimpleMMFFMinimization(**kwargs)
    
    def create_metal_aware_scoring_function(self, base_scoring_function: Any, 
                                          metal_centers: Optional[List] = None,
                                          metal_constraints: Optional[List] = None,
                                          metal_weight: float = 10.0):
        """Create metal-aware scoring function."""
        try:
            from .metal_docking import MetalDockingScorer
            return MetalDockingScorer(
                base_scoring_function=base_scoring_function,
                metal_centers=metal_centers or [],
                metal_constraints=metal_constraints or [],
                metal_weight=metal_weight
            )
        except ImportError:
            self.logger.warning("Metal scoring not available, using base scoring function")
            return base_scoring_function
    
    def get_available_algorithms(self) -> Dict[str, str]:
        """Get list of available algorithms with descriptions."""
        algorithms = {
            'genetic': 'Advanced genetic algorithm with tournament selection and adaptive mutation',
            'random': 'Random search with adaptive radius and local optimization', 
            'monte-carlo': 'Monte Carlo sampling with Metropolis criterion and simulated annealing',
            'pandadock': 'PANDADOCK algorithm with MD conformer generation and simulated annealing',
            'metal': 'Metal-aware docking with coordination constraints and specialized scoring'
        }
        return algorithms