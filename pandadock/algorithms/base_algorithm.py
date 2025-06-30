"""
Base algorithm interface for all docking algorithms.

This module defines the standard interface that all docking algorithms
must implement, ensuring consistency and interchangeability.
"""

from abc import ABC, abstractmethod
from typing import List, Tuple, Any, Optional, Dict
import logging


class BaseAlgorithm(ABC):
    """
    Abstract base class for all docking algorithms.
    
    All docking algorithms must inherit from this class and implement
    the required methods.
    """
    
    def __init__(self, scoring_function: Any, **kwargs):
        """
        Initialize the algorithm.
        
        Args:
            scoring_function: Scoring function to evaluate poses
            **kwargs: Algorithm-specific parameters
        """
        self.scoring_function = scoring_function
        self.logger = logging.getLogger(self.__class__.__name__)
        
        # Common parameters
        self.max_iterations = kwargs.get('max_iterations', 100)
        self.output_dir = kwargs.get('output_dir', None)
        
        # Results storage
        self.results = []
        self.iteration_count = 0
    
    @abstractmethod
    def search(self, protein: Any, ligand: Any) -> List[Tuple[Any, float]]:
        """
        Run the docking search algorithm.
        
        Args:
            protein: Protein molecule object
            ligand: Ligand molecule object
            
        Returns:
            List of (pose, score) tuples
        """
        pass
    
    def validate_inputs(self, protein: Any, ligand: Any) -> bool:
        """
        Validate input molecules.
        
        Args:
            protein: Protein molecule
            ligand: Ligand molecule
            
        Returns:
            True if inputs are valid
        """
        if protein is None:
            self.logger.error("Protein is None")
            return False
        
        if ligand is None:
            self.logger.error("Ligand is None")
            return False
        
        if not hasattr(protein, 'coords'):
            self.logger.error("Protein missing coordinates")
            return False
        
        if not hasattr(ligand, 'coords'):
            self.logger.error("Ligand missing coordinates")
            return False
        
        return True
    
    def get_algorithm_info(self) -> Dict[str, Any]:
        """
        Get information about the algorithm.
        
        Returns:
            Dictionary with algorithm information
        """
        return {
            'name': self.__class__.__name__,
            'max_iterations': self.max_iterations,
            'current_iteration': self.iteration_count,
            'total_results': len(self.results)
        }
    
    def cleanup(self) -> None:
        """Clean up algorithm resources."""
        pass