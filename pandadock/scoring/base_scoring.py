"""
Base scoring function interface.

This module defines the standard interface that all scoring functions
must implement.
"""

from abc import ABC, abstractmethod
from typing import List, Any, Dict
import logging


class BaseScoringFunction(ABC):
    """
    Abstract base class for all scoring functions.
    
    All scoring functions must inherit from this class and implement
    the required methods.
    """
    
    def __init__(self):
        """Initialize the scoring function."""
        self.logger = logging.getLogger(self.__class__.__name__)
        self.evaluation_count = 0
    
    @abstractmethod
    def score(self, protein: Any, ligand_pose: Any) -> float:
        """
        Score a single protein-ligand pose.
        
        Args:
            protein: Protein molecule object
            ligand_pose: Ligand pose to score
            
        Returns:
            Score value (lower is better)
        """
        pass
    
    def score_batch(self, protein: Any, ligand_poses: List[Any]) -> List[float]:
        """
        Score multiple poses efficiently.
        
        Args:
            protein: Protein molecule object
            ligand_poses: List of ligand poses to score
            
        Returns:
            List of score values
        """
        scores = []
        for pose in ligand_poses:
            score = self.score(protein, pose)
            scores.append(score)
            self.evaluation_count += 1
        
        return scores
    
    def score_with_components(self, protein: Any, ligand_pose: Any) -> Dict[str, float]:
        """
        Score a pose and return both total score and energy components.
        
        Args:
            protein: Protein molecule object
            ligand_pose: Ligand pose to score
            
        Returns:
            Dictionary with 'total_score' and energy components
        """
        total_score = self.score(protein, ligand_pose)
        return {
            'total_score': total_score,
            'van_der_waals': 0.0,
            'hydrogen_bonds': 0.0,
            'electrostatic': 0.0,
            'desolvation': 0.0,
            'hydrophobic': 0.0,
            'entropy': 0.0,
            'clash': 0.0
        }
    
    def get_scoring_info(self) -> Dict[str, Any]:
        """
        Get information about the scoring function.
        
        Returns:
            Dictionary with scoring function information
        """
        return {
            'name': self.__class__.__name__,
            'evaluation_count': self.evaluation_count
        }
    
    def reset_counters(self) -> None:
        """Reset evaluation counters."""
        self.evaluation_count = 0