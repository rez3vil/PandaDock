"""
Monte Carlo sampling algorithm for molecular docking.

This module implements Metropolis Monte Carlo sampling with optional
simulated annealing for pose optimization.
"""

import numpy as np
import logging
from typing import List, Tuple, Any, Optional, Dict
from copy import deepcopy
import time

from .base_algorithm import BaseAlgorithm
from ..utils_new.math_utils import rotation_matrix_from_euler, apply_rotation_translation

logger = logging.getLogger(__name__)


class MonteCarloAlgorithm(BaseAlgorithm):
    """Monte Carlo sampling algorithm with Metropolis criterion."""
    
    def __init__(self, scoring_function: Any, 
                 temperature: float = 300.0, n_steps: int = 1000,
                 max_translation: float = 2.0, max_rotation: float = 0.3,
                 cooling_factor: float = 0.95, annealing: bool = True,
                 equilibration_steps: int = 100, **kwargs):
        """
        Initialize Monte Carlo algorithm.
        
        Args:
            scoring_function: Scoring function to evaluate poses
            temperature: Simulation temperature (K)
            n_steps: Number of Monte Carlo steps
            max_translation: Maximum translation step size (Å)
            max_rotation: Maximum rotation step size (radians)
            cooling_factor: Temperature reduction factor for annealing
            annealing: Whether to use simulated annealing
            equilibration_steps: Number of equilibration steps
        """
        super().__init__(scoring_function, **kwargs)
        
        self.temperature = temperature
        self.n_steps = n_steps
        self.max_translation = max_translation
        self.max_rotation = max_rotation
        self.cooling_factor = cooling_factor
        self.annealing = annealing
        self.equilibration_steps = equilibration_steps
        
        self.logger = logging.getLogger(__name__)
        
        # Physical constants
        self.kb = 0.001987  # Boltzmann constant in kcal/(mol·K)
        
        # Statistics tracking
        self.accepted_moves = 0
        self.total_moves = 0
        self.energy_history = []
        
    def search(self, protein: Any, ligand: Any) -> List[Tuple[Any, float]]:
        """
        Perform Monte Carlo sampling.
        
        Args:
            protein: Protein object
            ligand: Ligand object
            
        Returns:
            List of (pose, score) tuples from the sampling trajectory
        """
        if not self.validate_inputs(protein, ligand):
            raise ValueError("Invalid protein or ligand inputs")
        
        self.logger.info(f"Starting Monte Carlo sampling")
        self.logger.info(f"Parameters: T={self.temperature}K, {self.n_steps} steps")
        self.logger.info(f"Step sizes: {self.max_translation}Å translation, {self.max_rotation:.2f} rad rotation")
        
        start_time = time.time()
        
        # Initialize starting pose
        current_pose = self._initialize_pose(protein, ligand)
        current_score = self.scoring_function.score(protein, current_pose)
        
        best_pose = deepcopy(current_pose)
        best_score = current_score
        
        # Reset statistics
        self.accepted_moves = 0
        self.total_moves = 0
        self.energy_history = []
        
        # Current temperature (for annealing)
        current_temp = self.temperature
        
        # Storage for trajectory
        trajectory = []
        save_frequency = max(1, self.n_steps // 100)  # Save ~100 poses
        
        self.logger.info("Starting Monte Carlo trajectory...")
        
        for step in range(self.n_steps):
            try:
                # Generate new pose
                new_pose = deepcopy(current_pose)
                self._apply_mc_move(new_pose)
                
                # Evaluate new pose
                new_score = self.scoring_function.score(protein, new_pose)
                
                # Metropolis acceptance criterion
                accept = self._metropolis_accept(current_score, new_score, current_temp)
                
                self.total_moves += 1
                
                if accept:
                    current_pose = new_pose
                    current_score = new_score
                    self.accepted_moves += 1
                    
                    # Update best if better
                    if new_score < best_score:
                        best_pose = deepcopy(new_pose)
                        best_score = new_score
                
                # Record energy
                self.energy_history.append(current_score)
                
                # Save pose to trajectory
                if step % save_frequency == 0 and step >= self.equilibration_steps:
                    trajectory.append((deepcopy(current_pose), current_score))
                
                # Temperature annealing
                if self.annealing and step % 50 == 0:
                    current_temp = max(self.temperature * 0.1, current_temp * self.cooling_factor)
                
                # Progress reporting
                if (step + 1) % 200 == 0:
                    acceptance_rate = self.accepted_moves / self.total_moves
                    self.logger.info(f"Step {step+1}/{self.n_steps}: "
                                   f"Current={current_score:.3f}, Best={best_score:.3f}, "
                                   f"Accept={acceptance_rate:.2%}, T={current_temp:.1f}K")
                
            except Exception as e:
                self.logger.warning(f"Error in MC step {step}: {e}")
                continue
        
        elapsed = time.time() - start_time
        final_acceptance = self.accepted_moves / max(self.total_moves, 1)
        
        self.logger.info(f"Monte Carlo completed in {elapsed:.1f}s")
        self.logger.info(f"Final statistics: {final_acceptance:.2%} acceptance rate")
        self.logger.info(f"Best score: {best_score:.3f}")
        
        # Add best pose to trajectory if not already there
        trajectory.append((best_pose, best_score))
        
        # Sort trajectory by score
        trajectory.sort(key=lambda x: x[1])
        
        return trajectory
    
    def _initialize_pose(self, protein: Any, ligand: Any) -> Any:
        """Initialize starting pose for Monte Carlo."""
        try:
            # Start with a random pose in the active site
            active_site = self._get_active_site(protein)
            if not active_site:
                return deepcopy(ligand)
            
            initial_pose = deepcopy(ligand)
            center = active_site['center']
            radius = active_site.get('radius', 10.0)
            
            # Random position within active site
            while True:
                random_point = np.random.uniform(-1, 1, 3)
                if np.linalg.norm(random_point) <= 1.0:
                    break
            
            position = center + random_point * radius * 0.5  # Start closer to center
            
            # Random orientation
            euler_angles = np.random.uniform(0, 2*np.pi, 3)
            rotation_matrix = rotation_matrix_from_euler(euler_angles)
            
            # Apply transformation
            self._apply_transformation(initial_pose, rotation_matrix, position)
            
            return initial_pose
            
        except Exception as e:
            self.logger.warning(f"Error initializing pose: {e}")
            return deepcopy(ligand)
    
    def _apply_mc_move(self, pose: Any) -> None:
        """Apply Monte Carlo move to pose."""
        try:
            # Random translation
            translation = np.random.uniform(-self.max_translation, self.max_translation, 3)
            
            # Random rotation around current centroid
            rotation_angles = np.random.uniform(-self.max_rotation, self.max_rotation, 3)
            rotation_matrix = rotation_matrix_from_euler(rotation_angles)
            
            # Get current centroid for rotation center
            centroid = self._get_pose_centroid(pose)
            
            if centroid is not None:
                # Rotate around centroid
                self._rotate_around_point(pose, rotation_matrix, centroid)
            
            # Apply translation
            self._translate_pose(pose, translation)
            
        except Exception as e:
            self.logger.debug(f"Error applying MC move: {e}")
    
    def _metropolis_accept(self, current_energy: float, new_energy: float, 
                          temperature: float) -> bool:
        """Metropolis acceptance criterion."""
        try:
            delta_e = new_energy - current_energy
            
            # Always accept better poses
            if delta_e <= 0:
                return True
            
            # Accept worse poses with Boltzmann probability
            boltzmann_factor = np.exp(-delta_e / (self.kb * temperature))
            return np.random.random() < boltzmann_factor
            
        except Exception as e:
            self.logger.debug(f"Error in Metropolis criterion: {e}")
            return False  # Reject on error
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get Monte Carlo statistics."""
        acceptance_rate = self.accepted_moves / max(self.total_moves, 1)
        
        stats = {
            'total_moves': self.total_moves,
            'accepted_moves': self.accepted_moves,
            'acceptance_rate': acceptance_rate,
            'n_energy_points': len(self.energy_history)
        }
        
        if self.energy_history:
            stats.update({
                'initial_energy': self.energy_history[0],
                'final_energy': self.energy_history[-1],
                'best_energy': min(self.energy_history),
                'average_energy': np.mean(self.energy_history),
                'energy_std': np.std(self.energy_history)
            })
        
        return stats
    
    def get_energy_trajectory(self) -> List[float]:
        """Get energy trajectory from sampling."""
        return self.energy_history.copy()
    
    # Helper methods for coordinate manipulation
    def _get_coordinates(self, molecule: Any) -> Optional[np.ndarray]:
        """Extract coordinates from molecule."""
        try:
            if hasattr(molecule, 'coords'):
                coords = np.array(molecule.coords)
            elif hasattr(molecule, 'coordinates'):
                coords = np.array(molecule.coordinates)
            elif hasattr(molecule, 'xyz'):
                coords = np.array(molecule.xyz)
            elif hasattr(molecule, 'atoms'):
                coords_list = []
                for atom in molecule.atoms:
                    atom_coords = self._get_atom_coords(atom)
                    if atom_coords is not None:
                        coords_list.append(atom_coords)
                if coords_list:
                    coords = np.array(coords_list)
                else:
                    return None
            else:
                return None
            
            return coords
            
        except Exception as e:
            self.logger.debug(f"Error getting coordinates: {e}")
            return None
    
    def _set_coordinates(self, molecule: Any, coords: np.ndarray) -> None:
        """Set coordinates in molecule."""
        try:
            if hasattr(molecule, 'coords'):
                molecule.coords = coords
            elif hasattr(molecule, 'coordinates'):
                molecule.coordinates = coords
            elif hasattr(molecule, 'xyz'):
                molecule.xyz = coords
            elif hasattr(molecule, 'atoms'):
                for i, atom in enumerate(molecule.atoms):
                    if i < len(coords):
                        if isinstance(atom, dict):
                            atom['coords'] = coords[i]
                        elif hasattr(atom, 'coords'):
                            atom.coords = coords[i]
        except Exception as e:
            self.logger.debug(f"Error setting coordinates: {e}")
    
    def _apply_transformation(self, pose: Any, rotation: np.ndarray, 
                            translation: np.ndarray) -> None:
        """Apply rotation and translation to pose."""
        try:
            coords = self._get_coordinates(pose)
            if coords is not None:
                if coords.ndim == 1:
                    coords = coords.reshape(-1, 3)
                
                transformed_coords = apply_rotation_translation(coords, rotation, translation)
                self._set_coordinates(pose, transformed_coords)
        except Exception as e:
            self.logger.debug(f"Error applying transformation: {e}")
    
    def _get_pose_centroid(self, pose: Any) -> Optional[np.ndarray]:
        """Calculate centroid of pose."""
        try:
            coords = self._get_coordinates(pose)
            if coords is not None:
                if coords.ndim == 1:
                    coords = coords.reshape(-1, 3)
                return np.mean(coords, axis=0)
            return None
        except Exception as e:
            self.logger.debug(f"Error calculating centroid: {e}")
            return None
    
    def _rotate_around_point(self, pose: Any, rotation_matrix: np.ndarray, 
                           center: np.ndarray) -> None:
        """Rotate pose around a specific point."""
        try:
            coords = self._get_coordinates(pose)
            if coords is not None:
                if coords.ndim == 1:
                    coords = coords.reshape(-1, 3)
                
                # Translate to origin, rotate, translate back
                centered_coords = coords - center
                rotated_coords = np.dot(centered_coords, rotation_matrix.T)
                final_coords = rotated_coords + center
                
                self._set_coordinates(pose, final_coords)
        except Exception as e:
            self.logger.debug(f"Error in rotation around point: {e}")
    
    def _translate_pose(self, pose: Any, translation: np.ndarray) -> None:
        """Apply translation to pose."""
        try:
            coords = self._get_coordinates(pose)
            if coords is not None:
                if coords.ndim == 1:
                    coords = coords.reshape(-1, 3)
                
                translated_coords = coords + translation
                self._set_coordinates(pose, translated_coords)
        except Exception as e:
            self.logger.debug(f"Error in translation: {e}")
    
    def _get_active_site(self, protein: Any) -> Optional[Dict]:
        """Extract active site information from protein."""
        try:
            if hasattr(protein, 'active_site') and protein.active_site:
                return protein.active_site
            elif hasattr(protein, 'binding_site'):
                return protein.binding_site
            else:
                # Create default active site
                atoms = self._get_protein_atoms(protein)
                if atoms:
                    coords = []
                    for atom in atoms:
                        coord = self._get_atom_coords(atom)
                        if coord is not None:
                            coords.append(coord)
                    
                    if coords:
                        coords_array = np.array(coords)
                        center = np.mean(coords_array, axis=0)
                        max_dist = np.max(np.linalg.norm(coords_array - center, axis=1))
                        
                        return {
                            'center': center,
                            'radius': max_dist + 10.0
                        }
                
                return None
        except Exception as e:
            self.logger.debug(f"Error getting active site: {e}")
            return None
    
    def _get_protein_atoms(self, protein: Any) -> List[Any]:
        """Get protein atoms."""
        try:
            if hasattr(protein, 'active_site') and protein.active_site and 'atoms' in protein.active_site:
                return protein.active_site['atoms']
            elif hasattr(protein, 'atoms'):
                return protein.atoms
            else:
                return []
        except:
            return []
    
    def _get_atom_coords(self, atom: Any) -> Optional[np.ndarray]:
        """Extract coordinates from atom."""
        try:
            if isinstance(atom, dict):
                coords = atom.get('coords', atom.get('coordinates', atom.get('xyz')))
                if coords is not None:
                    return np.array(coords)
            elif hasattr(atom, 'coords'):
                return np.array(atom.coords)
            elif hasattr(atom, 'coordinates'):
                return np.array(atom.coordinates)
            elif hasattr(atom, 'xyz'):
                return np.array(atom.xyz)
            
            return None
        except:
            return None


class MonteCarloSampling:
    """
    Standalone Monte Carlo sampling class for compatibility with old code.
    
    This provides the same interface as the original physics.py implementation.
    """
    
    def __init__(self, scoring_function: Any, temperature: float = 300.0, 
                 n_steps: int = 1000, max_translation: float = 2.0,
                 max_rotation: float = 0.3, cooling_factor: float = 0.95):
        """Initialize Monte Carlo sampling."""
        self.mc_algorithm = MonteCarloAlgorithm(
            scoring_function=scoring_function,
            temperature=temperature,
            n_steps=n_steps,
            max_translation=max_translation,
            max_rotation=max_rotation,
            cooling_factor=cooling_factor,
            annealing=True
        )
    
    def run_sampling(self, protein: Any, ligand: Any) -> List[Tuple[Any, float]]:
        """Run Monte Carlo sampling (compatibility method)."""
        return self.mc_algorithm.search(protein, ligand)