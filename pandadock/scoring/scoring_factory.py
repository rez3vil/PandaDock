"""
Factory for creating scoring functions.

This module provides a centralized way to create different scoring
functions with hardware-appropriate implementations.
"""

from typing import Any, Dict, Optional
import logging

from ..hardware import ComputeBackend
from .base_scoring import BaseScoringFunction


class ScoringFunctionFactory:
    """Factory for creating scoring functions."""
    
    def __init__(self, compute_backend: ComputeBackend = None):
        """Initialize the scoring function factory."""
        self.compute_backend = compute_backend
        self.logger = logging.getLogger(__name__)
    
    def create_scoring_function(self, scoring_type: str = 'standard',
                              enhanced: bool = False,
                              physics_based: bool = False,
                              **kwargs) -> BaseScoringFunction:
        """
        Create a scoring function.
        
        Args:
            scoring_type: Type of scoring function
            enhanced: Whether to use enhanced scoring
            physics_based: Whether to use physics-based scoring
            **kwargs: Additional parameters
            
        Returns:
            Configured scoring function instance
        """
        self.logger.info(f"Creating {scoring_type} scoring function")
        
        try:
            if physics_based:
                return self._create_physics_scoring(**kwargs)
            elif enhanced:
                return self._create_enhanced_scoring(**kwargs)
            else:
                return self._create_standard_scoring(**kwargs)
                
        except Exception as e:
            self.logger.error(f"Failed to create scoring function: {e}")
            # Fallback to basic scoring
            self.logger.info("Falling back to basic scoring function")
            return self._create_basic_scoring(**kwargs)
    
    def _create_standard_scoring(self, **kwargs) -> BaseScoringFunction:
        """Create standard composite scoring function."""
        try:
            from ..unified_scoring import CompositeScoringFunction
            
            # Adapter to make old scoring work with new interface
            class StandardScoringAdapter(BaseScoringFunction):
                def __init__(self, **kwargs):
                    super().__init__()
                    self.old_scoring = CompositeScoringFunction()
                
                def score(self, protein, ligand_pose):
                    return self.old_scoring.score(protein, ligand_pose)
                
                def score_batch(self, protein, ligand_poses):
                    return [self.score(protein, pose) for pose in ligand_poses]
            
            return StandardScoringAdapter(**kwargs)
            
        except ImportError:
            return self._create_basic_scoring(**kwargs)
    
    def _create_enhanced_scoring(self, **kwargs) -> BaseScoringFunction:
        """Create enhanced scoring function."""
        try:
            from ..unified_scoring import EnhancedScoringFunction
            
            class EnhancedScoringAdapter(BaseScoringFunction):
                def __init__(self, **kwargs):
                    super().__init__()
                    self.old_scoring = EnhancedScoringFunction()
                
                def score(self, protein, ligand_pose):
                    return self.old_scoring.score(protein, ligand_pose)
                
                def score_batch(self, protein, ligand_poses):
                    return [self.score(protein, pose) for pose in ligand_poses]
            
            return EnhancedScoringAdapter(**kwargs)
            
        except ImportError:
            self.logger.warning("Enhanced scoring not available, using standard")
            return self._create_standard_scoring(**kwargs)
    
    def _create_physics_scoring(self, **kwargs) -> BaseScoringFunction:
        """Create physics-based scoring function."""
        try:
            from .physics_based_scoring import PhysicsBasedScoringFunction
            return PhysicsBasedScoringFunction(**kwargs)
        except ImportError:
            self.logger.warning("Physics scoring not available, using enhanced")
            return self._create_enhanced_scoring(**kwargs)
    
    def _create_basic_scoring(self, **kwargs) -> BaseScoringFunction:
        """Create basic fallback scoring function."""
        
        class BasicScoringFunction(BaseScoringFunction):
            """Simple distance-based scoring function."""
            
            def score(self, protein, ligand_pose):
                """Simple scoring based on minimum distance."""
                try:
                    import numpy as np
                    
                    protein_coords = protein.coords
                    ligand_coords = ligand_pose.coords if hasattr(ligand_pose, 'coords') else ligand_pose
                    
                    # Calculate minimum distance
                    distances = np.linalg.norm(
                        protein_coords[:, np.newaxis, :] - ligand_coords[np.newaxis, :, :],
                        axis=2
                    )
                    min_distance = np.min(distances)
                    
                    # Simple scoring: penalize very close and very far
                    if min_distance < 1.0:
                        return 100.0  # Severe clash penalty
                    elif min_distance > 15.0:
                        return 25.0   # Too far penalty (reduced from 50.0)
                    elif min_distance > 10.0:
                        return 10.0   # Moderately far penalty
                    else:
                        return -min_distance  # Better when closer (within reason)
                        
                except Exception:
                    return 0.0  # Neutral score on error
            
            def score_batch(self, protein, ligand_poses):
                return [self.score(protein, pose) for pose in ligand_poses]
        
        return BasicScoringFunction(**kwargs)
    
    def get_available_scoring_functions(self) -> Dict[str, str]:
        """Get list of available scoring functions with descriptions."""
        scoring_functions = {
            'standard': 'Standard composite scoring function',
            'enhanced': 'Enhanced scoring with additional energy terms',
            'physics': 'Physics-based scoring with molecular mechanics and solvation',
            'basic': 'Simple distance-based scoring (fallback)'
        }
        return scoring_functions