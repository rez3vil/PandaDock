"""
Robust pose generation system for molecular docking.

This module provides reliable pose generation with error handling,
validation, and retry mechanisms to handle failures gracefully.
"""

import numpy as np
import logging
from typing import List, Tuple, Optional, Dict, Any
from dataclasses import dataclass
import time

from ..hardware import ComputeBackend, PerformanceMonitor


@dataclass
class PoseGenerationConfig:
    """Configuration for pose generation."""
    max_attempts: int = 5
    validation_enabled: bool = True
    min_distance_threshold: float = 1.0  # Minimum distance between atoms
    max_distance_threshold: float = 50.0  # Maximum reasonable distance
    energy_threshold: Optional[float] = None  # Maximum acceptable energy
    timeout_seconds: float = 60.0  # Timeout for pose generation


class PoseGenerationError(Exception):
    """Exception raised during pose generation."""
    pass


class PoseValidator:
    """Validates generated poses for geometric and energetic reasonableness."""
    
    def __init__(self, config: PoseGenerationConfig):
        """Initialize pose validator."""
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def validate_pose(self, pose_coords: np.ndarray, protein_coords: np.ndarray,
                     scoring_function: Optional[Any] = None) -> Tuple[bool, str]:
        """
        Validate a generated pose.
        
        Args:
            pose_coords: Ligand pose coordinates
            protein_coords: Protein coordinates
            scoring_function: Optional scoring function for energy validation
            
        Returns:
            Tuple of (is_valid, reason)
        """
        try:
            # Check for NaN or infinite coordinates
            if not np.isfinite(pose_coords).all():
                return False, "Invalid coordinates (NaN or infinite values)"
            
            # Check internal distances
            if not self._check_internal_distances(pose_coords):
                return False, "Invalid internal ligand geometry"
            
            # Check protein-ligand distances
            if not self._check_protein_ligand_distances(pose_coords, protein_coords):
                return False, "Invalid protein-ligand distances"
            
            # Check energy if scoring function provided
            if scoring_function and self.config.energy_threshold:
                energy = scoring_function.score_pose(pose_coords, protein_coords)
                if energy > self.config.energy_threshold:
                    return False, f"Energy too high: {energy:.2f}"
            
            return True, "Valid pose"
            
        except Exception as e:
            return False, f"Validation error: {str(e)}"
    
    def _check_internal_distances(self, coords: np.ndarray) -> bool:
        """Check internal ligand distances are reasonable."""
        if len(coords) < 2:
            return True
        
        # Calculate pairwise distances
        distances = np.linalg.norm(
            coords[:, np.newaxis, :] - coords[np.newaxis, :, :], 
            axis=2
        )
        
        # Remove diagonal (self-distances)
        np.fill_diagonal(distances, np.inf)
        
        # Check minimum distances
        min_distance = np.min(distances)
        if min_distance < self.config.min_distance_threshold:
            return False
        
        # Check maximum distances (for reasonable molecular size)
        max_distance = np.max(distances)
        if max_distance > self.config.max_distance_threshold:
            return False
        
        return True
    
    def _check_protein_ligand_distances(self, ligand_coords: np.ndarray, 
                                      protein_coords: np.ndarray) -> bool:
        """Check protein-ligand distances are reasonable."""
        if len(protein_coords) == 0 or len(ligand_coords) == 0:
            return True
        
        # Calculate minimum distance between ligand and protein
        distances = np.linalg.norm(
            ligand_coords[:, np.newaxis, :] - protein_coords[np.newaxis, :, :],
            axis=2
        )
        
        min_distance = np.min(distances)
        
        # Too close (severe clash)
        if min_distance < 0.5:
            return False
        
        # Too far (not interacting)
        if min_distance > 20.0:
            return False
        
        return True


class PoseGenerator:
    """
    Robust pose generation system with error handling and validation.
    
    Generates molecular poses with automatic retry mechanisms and
    validation to ensure high-quality docking results.
    """
    
    def __init__(self, compute_backend: ComputeBackend = None, 
                 performance_monitor: PerformanceMonitor = None,
                 config: Optional[PoseGenerationConfig] = None):
        """
        Initialize pose generator.
        
        Args:
            compute_backend: Compute backend for calculations (optional)
            performance_monitor: Performance monitoring system (optional)
            config: Pose generation configuration (optional)
        """
        self.compute_backend = compute_backend
        self.performance_monitor = performance_monitor
        self.config = config or PoseGenerationConfig()
        self.validator = PoseValidator(self.config)
        self.logger = logging.getLogger(__name__)
        
        # Statistics
        self.generation_stats = {
            'total_attempts': 0,
            'successful_generations': 0,
            'validation_failures': 0,
            'timeout_failures': 0,
            'error_failures': 0
        }
    
    def generate_initial_poses(self, ligand: Any, protein: Any, 
                             n_poses: int = 100) -> List[np.ndarray]:
        """
        Generate initial poses for docking search.
        
        Args:
            ligand: Ligand molecule object
            protein: Protein molecule object
            n_poses: Number of poses to generate
            
        Returns:
            List of valid pose coordinates
        """
        self.logger.info(f"Generating {n_poses} initial poses")
        
        valid_poses = []
        binding_site = self._get_binding_site(protein)
        
        with self.performance_monitor.start_operation("pose_generation"):
            for i in range(n_poses):
                pose = self._generate_single_pose(ligand, binding_site, protein)
                if pose is not None:
                    valid_poses.append(pose)
                    
                # Log progress periodically
                if (i + 1) % 20 == 0:
                    success_rate = len(valid_poses) / (i + 1)
                    self.logger.info(f"Generated {len(valid_poses)}/{i + 1} valid poses "
                                   f"(success rate: {success_rate:.1%})")
        
        self.logger.info(f"Generated {len(valid_poses)} valid poses out of {n_poses} attempts")
        
        if len(valid_poses) == 0:
            raise PoseGenerationError("No valid poses could be generated")
        
        return valid_poses
    
    def _generate_single_pose(self, ligand: Any, binding_site: Dict[str, Any], 
                            protein: Any) -> Optional[np.ndarray]:
        """
        Generate a single valid pose with retry mechanism.
        
        Args:
            ligand: Ligand molecule
            binding_site: Binding site information
            protein: Protein molecule
            
        Returns:
            Valid pose coordinates or None if generation failed
        """
        for attempt in range(self.config.max_attempts):
            self.generation_stats['total_attempts'] += 1
            
            try:
                # Generate pose with timeout
                start_time = time.time()
                pose_coords = self._generate_pose_attempt(ligand, binding_site)
                
                # Check timeout
                if time.time() - start_time > self.config.timeout_seconds:
                    self.generation_stats['timeout_failures'] += 1
                    continue
                
                # Validate pose if enabled
                if self.config.validation_enabled:
                    is_valid, reason = self.validator.validate_pose(
                        pose_coords, protein.coords
                    )
                    
                    if not is_valid:
                        self.generation_stats['validation_failures'] += 1
                        self.logger.debug(f"Pose validation failed: {reason}")
                        continue
                
                # Success
                self.generation_stats['successful_generations'] += 1
                return pose_coords
                
            except Exception as e:
                self.generation_stats['error_failures'] += 1
                self.logger.debug(f"Pose generation attempt {attempt + 1} failed: {e}")
                
                # For critical errors, stop trying
                if isinstance(e, MemoryError):
                    self.logger.error("Memory error during pose generation")
                    break
        
        return None
    
    def _generate_pose_attempt(self, ligand: Any, binding_site: Dict[str, Any]) -> np.ndarray:
        """
        Single attempt to generate a pose.
        
        Args:
            ligand: Ligand molecule
            binding_site: Binding site information
            
        Returns:
            Generated pose coordinates
        """
        # Get ligand coordinates
        ligand_coords = ligand.coords.copy()
        
        # Generate random translation within binding site
        if binding_site.get('center') is not None and binding_site.get('radius') is not None:
            center = np.array(binding_site['center'])
            radius = binding_site['radius']
            
            # Random translation within sphere
            direction = np.random.randn(3)
            direction = direction / np.linalg.norm(direction)
            distance = np.random.uniform(0, radius * 0.8)  # Stay within 80% of radius
            translation = center + direction * distance
        else:
            # Random translation in reasonable range
            translation = np.random.uniform(-10, 10, 3)
        
        # Apply translation
        translated_coords = ligand_coords + translation
        
        # Generate random rotation
        rotated_coords = self._apply_random_rotation(translated_coords)
        
        # Apply small random perturbation to avoid identical poses
        perturbation = np.random.normal(0, 0.1, rotated_coords.shape)
        final_coords = rotated_coords + perturbation
        
        return final_coords
    
    def _apply_random_rotation(self, coords: np.ndarray) -> np.ndarray:
        """Apply random rotation to coordinates."""
        # Center coordinates
        center = np.mean(coords, axis=0)
        centered_coords = coords - center
        
        # Generate random rotation matrix
        rotation_matrix = self._random_rotation_matrix()
        
        # Apply rotation
        rotated_coords = np.dot(centered_coords, rotation_matrix.T)
        
        # Restore center
        return rotated_coords + center
    
    def _random_rotation_matrix(self) -> np.ndarray:
        """Generate a random 3D rotation matrix."""
        # Generate random quaternion
        q = np.random.randn(4)
        q = q / np.linalg.norm(q)
        
        # Convert to rotation matrix
        w, x, y, z = q
        rotation_matrix = np.array([
            [1 - 2*y**2 - 2*z**2, 2*x*y - 2*w*z, 2*x*z + 2*w*y],
            [2*x*y + 2*w*z, 1 - 2*x**2 - 2*z**2, 2*y*z - 2*w*x],
            [2*x*z - 2*w*y, 2*y*z + 2*w*x, 1 - 2*x**2 - 2*y**2]
        ])
        
        return rotation_matrix
    
    def _get_binding_site(self, protein: Any) -> Dict[str, Any]:
        """Get binding site information from protein."""
        binding_site = {}
        
        if hasattr(protein, 'active_site') and protein.active_site:
            binding_site = protein.active_site
        else:
            # Use protein center as fallback
            if hasattr(protein, 'coords'):
                center = np.mean(protein.coords, axis=0)
                binding_site = {
                    'center': center,
                    'radius': 10.0  # Default radius
                }
        
        return binding_site
    
    def generate_conformers(self, ligand: Any, n_conformers: int = 10) -> List[np.ndarray]:
        """
        Generate ligand conformers using robust methods.
        
        Args:
            ligand: Ligand molecule
            n_conformers: Number of conformers to generate
            
        Returns:
            List of conformer coordinates
        """
        self.logger.info(f"Generating {n_conformers} ligand conformers")
        
        conformers = []
        
        with self.performance_monitor.start_operation("conformer_generation"):
            # Include original conformation
            conformers.append(ligand.coords.copy())
            
            # Generate additional conformers
            for i in range(n_conformers - 1):
                conformer = self._generate_conformer(ligand)
                if conformer is not None:
                    conformers.append(conformer)
        
        self.logger.info(f"Generated {len(conformers)} valid conformers")
        return conformers
    
    def _generate_conformer(self, ligand: Any) -> Optional[np.ndarray]:
        """Generate a single conformer."""
        try:
            # Start with original coordinates
            coords = ligand.coords.copy()
            
            # Apply small random perturbations to torsion angles
            # This is a simplified approach - in practice, would use
            # proper conformational sampling methods
            
            # Apply random perturbations
            perturbation = np.random.normal(0, 0.5, coords.shape)
            perturbed_coords = coords + perturbation
            
            # Validate conformer
            if self.config.validation_enabled:
                is_valid, _ = self.validator.validate_pose(perturbed_coords, np.array([]))
                if not is_valid:
                    return None
            
            return perturbed_coords
            
        except Exception as e:
            self.logger.debug(f"Conformer generation failed: {e}")
            return None
    
    def optimize_pose(self, pose_coords: np.ndarray, protein: Any, 
                     scoring_function: Any) -> np.ndarray:
        """
        Optimize a pose using local optimization.
        
        Args:
            pose_coords: Initial pose coordinates
            protein: Protein molecule
            scoring_function: Scoring function for optimization
            
        Returns:
            Optimized pose coordinates
        """
        try:
            with self.performance_monitor.start_operation("pose_optimization"):
                # Use compute backend for optimization
                optimization_params = {
                    'scoring_function': scoring_function,
                    'protein_coords': protein.coords,
                    'max_iterations': 100,
                    'convergence_threshold': 0.01
                }
                
                optimized_coords = self.compute_backend.optimize_pose(
                    pose_coords, optimization_params
                )
                
                # Validate optimized pose
                if self.config.validation_enabled:
                    is_valid, reason = self.validator.validate_pose(
                        optimized_coords, protein.coords, scoring_function
                    )
                    
                    if not is_valid:
                        self.logger.warning(f"Optimized pose validation failed: {reason}")
                        return pose_coords  # Return original if optimization failed
                
                return optimized_coords
                
        except Exception as e:
            self.logger.warning(f"Pose optimization failed: {e}")
            return pose_coords  # Return original on failure
    
    def get_generation_statistics(self) -> Dict[str, Any]:
        """Get pose generation statistics."""
        total = self.generation_stats['total_attempts']
        if total == 0:
            return {'no_attempts': True}
        
        return {
            'total_attempts': total,
            'successful_generations': self.generation_stats['successful_generations'],
            'success_rate': self.generation_stats['successful_generations'] / total,
            'validation_failure_rate': self.generation_stats['validation_failures'] / total,
            'timeout_failure_rate': self.generation_stats['timeout_failures'] / total,
            'error_failure_rate': self.generation_stats['error_failures'] / total
        }
    
    def reset_statistics(self) -> None:
        """Reset generation statistics."""
        self.generation_stats = {
            'total_attempts': 0,
            'successful_generations': 0,
            'validation_failures': 0,
            'timeout_failures': 0,
            'error_failures': 0
        }