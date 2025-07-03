"""
Random search algorithm for molecular docking.

This module implements a robust random search algorithm with adaptive radius,
clash detection, and optional local optimization capabilities.
"""

import numpy as np
import logging
from typing import List, Tuple, Any, Optional, Dict, Union
from copy import deepcopy
import time
from tqdm import tqdm

from .base_algorithm import BaseAlgorithm
from ..utils_new.math_utils import rotation_matrix_from_euler, apply_rotation_translation

logger = logging.getLogger(__name__)


class RandomSearchAlgorithm(BaseAlgorithm):
    """Random search algorithm with adaptive radius and local optimization."""
    
    def __init__(self, scoring_function: Any, max_iterations: int = 1000,
                 initial_radius: float = 15.0, min_radius: float = 3.0,
                 radius_shrink_factor: float = 0.95, local_optimization: bool = False,
                 clash_threshold: float = 1.5, max_clashes: int = 5,
                 convergence_patience: int = 100, **kwargs):
        """
        Initialize random search algorithm.
        
        Args:
            scoring_function: Scoring function to evaluate poses
            max_iterations: Maximum number of search iterations
            initial_radius: Initial search radius around active site
            min_radius: Minimum search radius
            radius_shrink_factor: Factor to shrink radius when poses are found
            local_optimization: Whether to perform local optimization
            clash_threshold: Distance threshold for clash detection (1.5Å recommended)
            max_clashes: Maximum allowed clashes per pose
            convergence_patience: Iterations without improvement before stopping
        """
        super().__init__(scoring_function, **kwargs)
        
        self.max_iterations = max_iterations
        self.initial_radius = initial_radius
        self.min_radius = min_radius
        self.radius_shrink_factor = radius_shrink_factor
        self.local_optimization = local_optimization
        self.clash_threshold = clash_threshold
        self.max_clashes = max_clashes
        self.convergence_patience = convergence_patience
        
        self.logger = logging.getLogger(__name__)
        
        # Search state
        self.current_radius = initial_radius
        self.best_score = float('inf')
        self.iterations_without_improvement = 0
        self.pose_cache = []
        
    def search(self, protein: Any, ligand: Any) -> List[Tuple[Any, float]]:
        """
        Perform random search docking.
        
        Args:
            protein: Protein object
            ligand: Ligand object
            
        Returns:
            List of (pose, score) tuples sorted by score
        """
        if not self.validate_inputs(protein, ligand):
            raise ValueError("Invalid protein or ligand inputs")
        
        self.logger.info(f"Starting random search with {self.max_iterations} iterations")
        
        # Reset search state
        self.current_radius = self.initial_radius
        self.best_score = float('inf')
        self.iterations_without_improvement = 0
        self.pose_cache = []
        
        # Get active site information
        active_site = self._get_active_site(protein)
        if not active_site:
            raise ValueError("No active site defined for protein")
        
        center = active_site['center']
        max_radius = active_site.get('radius', self.initial_radius)
        self.current_radius = min(self.initial_radius, max_radius)
        
        start_time = time.time()
        results = []
        
        with tqdm(range(self.max_iterations), desc="Random Search", unit="iter") as pbar:
            for iteration in pbar:
                try:
                    # Generate random pose
                    pose = self._generate_random_pose(ligand, center, self.current_radius)
                    
                    # Check for clashes
                    if not self._is_pose_valid(protein, pose):
                        continue
                    
                    # Evaluate pose with energy components if available
                    if hasattr(self.scoring_function, 'score_with_components'):
                        score_data = self.scoring_function.score_with_components(protein, pose)
                        score = score_data['total_score']
                        pose.energy_components = {
                            'van_der_waals': score_data.get('van_der_waals', 0.0),
                            'hydrogen_bonds': score_data.get('hydrogen_bonds', 0.0),
                            'electrostatic': score_data.get('electrostatic', 0.0),
                            'desolvation': score_data.get('desolvation', 0.0),
                            'hydrophobic': score_data.get('hydrophobic', 0.0),
                            'entropy': score_data.get('entropy', 0.0),
                            'clash': score_data.get('clash', 0.0)
                        }
                    else:
                        # Fallback to simple scoring
                        score = self.scoring_function.score(protein, pose)
                    
                    # Store result
                    results.append((pose, score))
                    
                    # Update best score and adapt search
                    if score < self.best_score:
                        self.best_score = score
                        self.iterations_without_improvement = 0
                        
                        # Shrink radius when better poses are found
                        if self.current_radius > self.min_radius:
                            self.current_radius = max(
                                self.min_radius,
                                self.current_radius * self.radius_shrink_factor
                            )
                            self.logger.debug(f"Shrunk search radius to {self.current_radius:.2f}Å")
                    else:
                        self.iterations_without_improvement += 1
                    
                    # Update progress bar
                    pbar.set_postfix({
                        'best': f'{self.best_score:.3f}',
                        'poses': len(results),
                        'radius': f'{self.current_radius:.2f}Å'
                    })
                    
                    # Check convergence
                    if self.iterations_without_improvement >= self.convergence_patience:
                        self.logger.info(f"Converged after {iteration+1} iterations")
                        pbar.set_description("Random Search (Converged)")
                        pbar.close()
                        break
                    
                    # Progress reporting (less frequent to avoid cluttering with progress bar)
                    if (iteration + 1) % 200 == 0:
                        elapsed = time.time() - start_time
                        self.logger.info(f"Iteration {iteration+1}/{self.max_iterations}, "
                                       f"Best score: {self.best_score:.3f}, "
                                       f"Radius: {self.current_radius:.2f}Å, "
                                       f"Time: {elapsed:.1f}s")
                    
                except Exception as e:
                    self.logger.warning(f"Error in iteration {iteration}: {e}")
                    continue
        
        elapsed = time.time() - start_time
        self.logger.info(f"Random search completed: {len(results)} poses in {elapsed:.1f}s")
        
        # Sort results by score
        results.sort(key=lambda x: x[1])
        
        # Apply local optimization if requested
        if self.local_optimization and results:
            results = self._apply_local_optimization(protein, results)
        
        return results
    
    def _generate_random_pose(self, ligand: Any, center: np.ndarray, 
                            radius: float) -> Any:
        """Generate a random pose within the search radius."""
        # Create a copy of the ligand to modify
        pose = deepcopy(ligand)
        
        # Calculate ligand center of mass
        ligand_center = np.mean(pose.coords, axis=0)
        
        # Random translation within sphere
        # Use rejection sampling to ensure uniform distribution in sphere
        while True:
            random_point = np.random.uniform(-1, 1, 3)
            if np.linalg.norm(random_point) <= 1.0:
                break
        
        # Generate target position within sphere
        distance = radius * np.random.random()**(1/3)  # Uniform distribution in sphere
        target_position = center + random_point * distance
        
        # Calculate translation needed to move ligand center to target position
        translation = target_position - ligand_center
        
        # Random rotation
        euler_angles = np.random.uniform(0, 2*np.pi, 3)
        rotation_matrix = rotation_matrix_from_euler(euler_angles)
        
        # Apply transformation to pose
        if hasattr(pose, 'apply_transformation'):
            pose.apply_transformation(rotation_matrix, translation)
        else:
            # Fallback implementation
            self._apply_transformation_to_pose(pose, rotation_matrix, translation)
        
        return pose
    
    def _apply_transformation_to_pose(self, pose: Any, rotation: np.ndarray, 
                                    translation: np.ndarray) -> None:
        """Apply rotation and translation to pose coordinates."""
        try:
            # Try different coordinate attribute names
            coords = None
            coord_attr = None
            
            for attr in ['coords', 'coordinates', 'xyz', 'atoms']:
                if hasattr(pose, attr):
                    coords = getattr(pose, attr)
                    coord_attr = attr
                    break
            
            if coords is None:
                raise ValueError("Could not find coordinates in pose")
            
            # Handle different coordinate formats
            if coord_attr == 'atoms' and isinstance(coords, list):
                # Atom list format
                for atom in coords:
                    if isinstance(atom, dict) and 'coords' in atom:
                        atom_coords = np.array(atom['coords'])
                        atom_coords = apply_rotation_translation(
                            atom_coords.reshape(1, -1), rotation, translation
                        )[0]
                        atom['coords'] = atom_coords
                    elif hasattr(atom, 'coords'):
                        atom_coords = np.array(atom.coords)
                        atom_coords = apply_rotation_translation(
                            atom_coords.reshape(1, -1), rotation, translation
                        )[0]
                        atom.coords = atom_coords
            
            elif isinstance(coords, np.ndarray) and coords.ndim == 2:
                # Matrix format (N_atoms x 3)
                transformed_coords = apply_rotation_translation(coords, rotation, translation)
                setattr(pose, coord_attr, transformed_coords)
            
            elif isinstance(coords, np.ndarray) and coords.ndim == 1:
                # Flattened format
                coords_3d = coords.reshape(-1, 3)
                transformed_coords = apply_rotation_translation(coords_3d, rotation, translation)
                setattr(pose, coord_attr, transformed_coords.flatten())
            
            else:
                # Try to convert to array
                coords_array = np.array(coords)
                if coords_array.size > 0:
                    coords_3d = coords_array.reshape(-1, 3)
                    transformed_coords = apply_rotation_translation(coords_3d, rotation, translation)
                    setattr(pose, coord_attr, transformed_coords)
                else:
                    raise ValueError("Empty or invalid coordinates")
            
        except Exception as e:
            self.logger.warning(f"Failed to apply transformation to pose: {e}")
            # Don't raise - just continue with original pose
    
    def _is_pose_valid(self, protein: Any, pose: Any) -> bool:
        """Check if pose is valid (no severe clashes)."""
        try:
            # Get protein atoms
            protein_atoms = self._get_protein_atoms(protein)
            pose_atoms = self._get_pose_atoms(pose)
            
            if not protein_atoms or not pose_atoms:
                return True  # No atoms to check
            
            clash_count = 0
            
            # Check distances between protein and ligand atoms
            for p_atom in protein_atoms:
                p_coords = self._get_atom_coords(p_atom)
                if p_coords is None:
                    continue
                
                for l_atom in pose_atoms:
                    l_coords = self._get_atom_coords(l_atom)
                    if l_coords is None:
                        continue
                    
                    distance = np.linalg.norm(p_coords - l_coords)
                    
                    if distance < self.clash_threshold:
                        clash_count += 1
                        if clash_count > self.max_clashes:
                            return False
            
            return True
            
        except Exception as e:
            self.logger.debug(f"Error in pose validation: {e}")
            return True  # Default to valid if check fails
    
    def _get_protein_atoms(self, protein: Any) -> List[Any]:
        """Extract protein atoms for clash checking."""
        try:
            # First try to get atoms from active site
            if hasattr(protein, 'active_site') and protein.active_site and 'atoms' in protein.active_site:
                return protein.active_site['atoms']
            # Try atoms attribute
            elif hasattr(protein, 'atoms'):
                return protein.atoms
            # For ProteinStructure objects, create atoms list from coords and residues
            elif hasattr(protein, 'coords') and hasattr(protein, 'residues'):
                atoms = []
                for res_id, res_atoms in protein.residues.items():
                    for atom in res_atoms:
                        atoms.append(atom)
                return atoms
            # Fallback to creating atom-like objects from coords
            elif hasattr(protein, 'coords'):
                coords = protein.coords
                atoms = []
                for i, coord in enumerate(coords):
                    atom = {'coords': coord, 'atom_index': i}
                    atoms.append(atom)
                return atoms
            else:
                return []
        except Exception as e:
            self.logger.debug(f"Error getting protein atoms: {e}")
            return []
    
    def _get_pose_atoms(self, pose: Any) -> List[Any]:
        """Extract pose atoms for clash checking."""
        try:
            # Try atoms attribute first
            if hasattr(pose, 'atoms'):
                return pose.atoms
            elif hasattr(pose, 'molecule') and hasattr(pose.molecule, 'atoms'):
                return pose.molecule.atoms
            # For ligand objects with coords, create atom-like objects
            elif hasattr(pose, 'coords'):
                coords = pose.coords
                atoms = []
                for i, coord in enumerate(coords):
                    atom = {'coords': coord, 'atom_index': i}
                    atoms.append(atom)
                return atoms
            else:
                return []
        except Exception as e:
            self.logger.debug(f"Error getting pose atoms: {e}")
            return []
    
    def _get_atom_coords(self, atom: Any) -> Optional[np.ndarray]:
        """Extract coordinates from atom object."""
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
    
    def _get_active_site(self, protein: Any) -> Optional[Dict]:
        """Extract active site information from protein."""
        try:
            # Check for existing active site
            if hasattr(protein, 'active_site') and protein.active_site:
                return protein.active_site
            elif hasattr(protein, 'binding_site'):
                return protein.binding_site
            else:
                # Create default active site from protein coordinates
                coords = None
                
                # Try to get coordinates directly from protein
                if hasattr(protein, 'coords'):
                    coords = protein.coords
                elif hasattr(protein, 'get_center_of_mass'):
                    # Use center of mass method if available
                    center = protein.get_center_of_mass()
                    return {
                        'center': center,
                        'radius': 15.0  # Default radius
                    }
                else:
                    # Fallback to extracting from atoms
                    atoms = self._get_protein_atoms(protein)
                    if atoms:
                        atom_coords = []
                        for atom in atoms:
                            coord = self._get_atom_coords(atom)
                            if coord is not None:
                                atom_coords.append(coord)
                        
                        if atom_coords:
                            coords = np.array(atom_coords)
                
                # Calculate center and radius from coordinates
                if coords is not None and len(coords) > 0:
                    coords = np.array(coords)
                    center = np.mean(coords, axis=0)
                    max_dist = np.max(np.linalg.norm(coords - center, axis=1))
                    
                    return {
                        'center': center,
                        'radius': max_dist + 5.0  # Add some padding
                    }
                
                return None
        except Exception as e:
            self.logger.debug(f"Error getting active site: {e}")
            return None
    
    def _apply_local_optimization(self, protein: Any, 
                                results: List[Tuple[Any, float]]) -> List[Tuple[Any, float]]:
        """Apply local optimization to top poses."""
        try:
            self.logger.info("Applying local optimization to top poses")
            
            optimized_results = []
            n_optimize = min(5, len(results))  # Optimize top 5 poses (reduced from 10)
            
            for i, (pose, score) in enumerate(results[:n_optimize]):
                try:
                    optimized_pose, optimized_score = self._local_minimize(protein, pose, score)
                    optimized_results.append((optimized_pose, optimized_score))
                    
                    if optimized_score < score:
                        improvement = score - optimized_score
                        self.logger.debug(f"Pose {i+1}: improved by {improvement:.3f}")
                    
                except Exception as e:
                    self.logger.warning(f"Local optimization failed for pose {i+1}: {e}")
                    optimized_results.append((pose, score))
            
            # Add remaining poses without optimization
            optimized_results.extend(results[n_optimize:])
            
            # Re-sort by score
            optimized_results.sort(key=lambda x: x[1])
            
            return optimized_results
            
        except Exception as e:
            self.logger.error(f"Error in local optimization: {e}")
            return results
    
    def _local_minimize(self, protein: Any, pose: Any, 
                       initial_score: float) -> Tuple[Any, float]:
        """Perform local minimization of a pose."""
        try:
            best_pose = deepcopy(pose)
            best_score = initial_score
            
            # Simple gradient-free optimization using random perturbations
            n_steps = 10  # Reduced from 20 to 10 steps
            step_size = 0.5  # Angstroms
            angle_step = 0.1  # radians
            steps_without_improvement = 0
            max_no_improvement = 3  # Early stopping after 3 steps without improvement
            
            for step in range(n_steps):
                # Generate small random perturbation
                translation_delta = np.random.normal(0, step_size, 3)
                rotation_delta = np.random.normal(0, angle_step, 3)
                
                # Apply perturbation
                test_pose = deepcopy(best_pose)
                rotation_matrix = rotation_matrix_from_euler(rotation_delta)
                
                # Get current centroid for rotation center
                centroid = self._get_pose_centroid(test_pose)
                if centroid is not None:
                    self._rotate_around_point(test_pose, rotation_matrix, centroid)
                    self._translate_pose(test_pose, translation_delta)
                
                # Check if still valid
                if not self._is_pose_valid(protein, test_pose):
                    continue
                
                # Evaluate new pose
                new_score = self.scoring_function.score(protein, test_pose)
                
                # Accept if better
                if new_score < best_score:
                    best_pose = test_pose
                    best_score = new_score
                    step_size *= 1.1  # Increase step size on success
                    steps_without_improvement = 0  # Reset counter
                else:
                    step_size *= 0.9  # Decrease step size on failure
                    steps_without_improvement += 1
                
                # Early stopping if no improvement
                if steps_without_improvement >= max_no_improvement:
                    self.logger.debug(f"Early stopping local optimization at step {step+1}")
                    break
                
                # Ensure step size doesn't get too small or large
                step_size = max(0.1, min(2.0, step_size))
            
            return best_pose, best_score
            
        except Exception as e:
            self.logger.debug(f"Local minimization failed: {e}")
            return pose, initial_score
    
    def _get_pose_centroid(self, pose: Any) -> Optional[np.ndarray]:
        """Calculate centroid of pose coordinates."""
        try:
            atoms = self._get_pose_atoms(pose)
            if not atoms:
                return None
            
            coords = []
            for atom in atoms:
                coord = self._get_atom_coords(atom)
                if coord is not None:
                    coords.append(coord)
            
            if coords:
                return np.mean(coords, axis=0)
            else:
                return None
        except:
            return None
    
    def _rotate_around_point(self, pose: Any, rotation_matrix: np.ndarray, 
                           center: np.ndarray) -> None:
        """Rotate pose around a specific point."""
        try:
            # Translate to origin, rotate, translate back
            self._translate_pose(pose, -center)
            self._apply_transformation_to_pose(pose, rotation_matrix, np.zeros(3))
            self._translate_pose(pose, center)
        except Exception as e:
            self.logger.debug(f"Rotation around point failed: {e}")
    
    def _translate_pose(self, pose: Any, translation: np.ndarray) -> None:
        """Apply translation to pose."""
        try:
            self._apply_transformation_to_pose(pose, np.eye(3), translation)
        except Exception as e:
            self.logger.debug(f"Translation failed: {e}")


class SimpleRandomSearch(BaseAlgorithm):
    """Simplified random search for fallback compatibility."""
    
    def __init__(self, scoring_function: Any, max_iterations: int = 100, **kwargs):
        super().__init__(scoring_function, **kwargs)
        self.max_iterations = max_iterations
        self.logger = logging.getLogger(__name__)
    
    def search(self, protein: Any, ligand: Any) -> List[Tuple[Any, float]]:
        """Simple random search implementation."""
        if not self.validate_inputs(protein, ligand):
            raise ValueError("Invalid inputs")
        
        self.logger.info(f"Running simple random search with {self.max_iterations} iterations")
        
        results = []
        
        with tqdm(range(self.max_iterations), desc="Simple Random Search", unit="iter") as pbar:
            for i in pbar:
                try:
                    # Create random pose (very basic implementation)
                    pose = deepcopy(ligand)
                    
                    # Random translation with larger spread
                    translation = np.random.normal(0, 10.0, 3)
                    
                    # Random rotation
                    euler = np.random.uniform(0, 2*np.pi, 3)
                    rotation = rotation_matrix_from_euler(euler)
                    
                    # Apply transformation (basic attempt)
                    if hasattr(pose, 'coords'):
                        coords = np.array(pose.coords)
                        if coords.ndim == 2:
                            transformed = apply_rotation_translation(coords, rotation, translation)
                            pose.coords = transformed
                    
                    # Basic clash check before scoring
                    if self._basic_clash_check(protein, pose):
                        continue
                    
                    # Score pose
                    score = self.scoring_function.score(protein, pose)
                    results.append((pose, score))
                    
                    # Update progress bar
                    if results:
                        best_score = min(result[1] for result in results)
                        pbar.set_postfix({'best': f'{best_score:.3f}', 'poses': len(results)})
                    
                except Exception as e:
                    self.logger.debug(f"Error in simple random search iteration {i}: {e}")
                    continue
        
        # Sort by score
        results.sort(key=lambda x: x[1])
        
        self.logger.info(f"Simple random search completed: {len(results)} poses")
        return results
    
    def _basic_clash_check(self, protein: Any, pose: Any) -> bool:
        """Basic clash check for SimpleRandomSearch."""
        try:
            # Get coordinates using the same methods as the main class
            protein_coords = None
            pose_coords = None
            
            # Extract protein coordinates
            if hasattr(protein, 'coords'):
                protein_coords = protein.coords
            else:
                return False  # No clash if no coordinates
            
            # Extract pose coordinates  
            if hasattr(pose, 'coords'):
                pose_coords = pose.coords
            else:
                return False  # No clash if no coordinates
            
            if protein_coords is None or pose_coords is None:
                return False
            
            # Convert to numpy arrays
            protein_coords = np.array(protein_coords)
            pose_coords = np.array(pose_coords)
            
            # Check minimum distance
            if protein_coords.size > 0 and pose_coords.size > 0:
                # Ensure proper shape (N, 3)
                if protein_coords.ndim == 1:
                    protein_coords = protein_coords.reshape(-1, 3)
                if pose_coords.ndim == 1:
                    pose_coords = pose_coords.reshape(-1, 3)
                
                # Calculate all pairwise distances
                distances = np.linalg.norm(
                    protein_coords[:, np.newaxis, :] - pose_coords[np.newaxis, :, :],
                    axis=2
                )
                min_distance = np.min(distances)
                
                # Return True if clash detected (distance < 1.5Å)
                return min_distance < 1.5
            
            return False
        except Exception as e:
            self.logger.debug(f"Clash check failed: {e}")
            return False  # Default to no clash if check fails