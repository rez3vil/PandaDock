"""
PANDADOCK algorithm for molecular docking.

This module implements the PANDADOCK algorithm which uses high-temperature
molecular dynamics followed by simulated annealing for pose generation.
"""

import numpy as np
import logging
from typing import List, Tuple, Any, Optional, Dict
from copy import deepcopy
import time

from .base_algorithm import BaseAlgorithm
from ..utils_new.math_utils import rotation_matrix_from_euler, apply_rotation_translation

logger = logging.getLogger(__name__)


class PANDADOCKAlgorithm(BaseAlgorithm):
    """PANDADOCK algorithm with MD conformer generation and simulated annealing."""
    
    def __init__(self, scoring_function: Any, 
                 high_temp: float = 1000.0, target_temp: float = 300.0,
                 num_conformers: int = 10, num_orientations: int = 10,
                 md_steps: int = 1000, minimize_steps: int = 200,
                 cooling_factor: float = 0.95, **kwargs):
        """
        Initialize PANDADOCK algorithm.
        
        Args:
            scoring_function: Scoring function to evaluate poses
            high_temp: High temperature for MD conformer generation (K)
            target_temp: Final target temperature (K)
            num_conformers: Number of ligand conformers to generate
            num_orientations: Number of orientations per conformer
            md_steps: Number of simulated annealing steps
            minimize_steps: Number of final minimization steps
            cooling_factor: Temperature reduction factor per cycle
        """
        super().__init__(scoring_function, **kwargs)
        
        self.high_temp = high_temp
        self.target_temp = target_temp
        self.num_conformers = num_conformers
        self.num_orientations = num_orientations
        self.md_steps = md_steps
        self.minimize_steps = minimize_steps
        self.cooling_factor = cooling_factor
        
        self.logger = logging.getLogger(__name__)
        
        # Physical constants
        self.kb = 0.001987  # Boltzmann constant in kcal/(mol·K)
        
    def search(self, protein: Any, ligand: Any) -> List[Tuple[Any, float]]:
        """
        Perform PANDADOCK docking protocol.
        
        Args:
            protein: Protein object
            ligand: Ligand object
            
        Returns:
            List of (pose, score) tuples sorted by score
        """
        if not self.validate_inputs(protein, ligand):
            raise ValueError("Invalid protein or ligand inputs")
        
        self.logger.info("Starting PANDADOCK algorithm")
        self.logger.info(f"Protocol: {self.num_conformers} conformers × {self.num_orientations} orientations")
        self.logger.info(f"Temperature: {self.high_temp}K → {self.target_temp}K")
        
        start_time = time.time()
        all_results = []
        
        # Get active site information
        active_site = self._get_active_site(protein)
        if not active_site:
            raise ValueError("No active site defined for protein")
        
        # Step 1: Generate ligand conformers using high-temperature MD
        self.logger.info("Step 1: Generating ligand conformers...")
        conformers = self._generate_conformers(ligand)
        self.logger.info(f"Generated {len(conformers)} conformers")
        
        # Step 2: For each conformer, generate random orientations and positions
        for conf_idx, conformer in enumerate(conformers):
            self.logger.info(f"Processing conformer {conf_idx + 1}/{len(conformers)}")
            
            orientations = self._generate_orientations(conformer, active_site, self.num_orientations)
            
            # Step 3: Simulated annealing for each orientation
            for orient_idx, oriented_pose in enumerate(orientations):
                try:
                    self.logger.debug(f"  Annealing orientation {orient_idx + 1}/{len(orientations)}")
                    
                    annealed_pose = self._simulated_annealing(protein, oriented_pose)
                    
                    # Step 4: Final minimization
                    minimized_pose = self._final_minimization(protein, annealed_pose)
                    
                    # Evaluate final pose
                    final_score = self.scoring_function.score(protein, minimized_pose)
                    
                    # Check pose validity
                    if self._check_pose_validity(protein, minimized_pose, final_score):
                        all_results.append((minimized_pose, final_score))
                    
                except Exception as e:
                    self.logger.warning(f"Error processing conformer {conf_idx}, orientation {orient_idx}: {e}")
                    continue
        
        elapsed = time.time() - start_time
        self.logger.info(f"PANDADOCK completed: {len(all_results)} poses in {elapsed:.1f}s")
        
        # Sort results by score
        all_results.sort(key=lambda x: x[1])
        
        return all_results
    
    def _generate_conformers(self, ligand: Any) -> List[Any]:
        """Generate ligand conformers using high-temperature MD simulation."""
        conformers = []
        
        try:
            # Start with original ligand
            conformers.append(deepcopy(ligand))
            
            # Generate additional conformers by perturbing the structure
            base_conformer = deepcopy(ligand)
            
            for i in range(self.num_conformers - 1):
                # Create perturbed conformer
                new_conformer = deepcopy(base_conformer)
                
                # Apply high-temperature MD-like perturbations
                self._apply_md_perturbation(new_conformer)
                
                conformers.append(new_conformer)
            
            return conformers
            
        except Exception as e:
            self.logger.warning(f"Error generating conformers: {e}")
            return [deepcopy(ligand)]  # Return at least the original
    
    def _apply_md_perturbation(self, conformer: Any) -> None:
        """Apply MD-like perturbations to generate new conformers."""
        try:
            # Get coordinates
            coords = self._get_coordinates(conformer)
            if coords is None or coords.size == 0:
                return
            
            # Ensure coords is 2D (n_atoms, 3)
            if coords.ndim == 1:
                coords = coords.reshape(-1, 3)
            
            # Calculate thermal energy at high temperature
            thermal_energy = self.kb * self.high_temp
            
            # Apply random perturbations based on thermal motion
            # Displacement scale based on thermal energy (simplified)
            displacement_scale = np.sqrt(thermal_energy) * 0.5  # Angstroms
            
            # Generate random displacements
            random_displacements = np.random.normal(0, displacement_scale, coords.shape)
            
            # Apply perturbations
            perturbed_coords = coords + random_displacements
            
            # Update conformer coordinates
            self._set_coordinates(conformer, perturbed_coords)
            
        except Exception as e:
            self.logger.debug(f"Error applying MD perturbation: {e}")
    
    def _generate_orientations(self, conformer: Any, active_site: Dict, 
                             num_orientations: int) -> List[Any]:
        """Generate random orientations and positions for a conformer."""
        orientations = []
        
        center = active_site['center']
        radius = active_site.get('radius', 10.0)
        
        for i in range(num_orientations):
            try:
                # Create copy of conformer
                oriented_pose = deepcopy(conformer)
                
                # Random position within active site sphere
                # Use rejection sampling for uniform distribution
                while True:
                    random_point = np.random.uniform(-1, 1, 3)
                    if np.linalg.norm(random_point) <= 1.0:
                        break
                
                position = center + random_point * radius * np.random.random()
                
                # Random orientation
                euler_angles = np.random.uniform(0, 2*np.pi, 3)
                rotation_matrix = rotation_matrix_from_euler(euler_angles)
                
                # Apply transformation
                self._apply_transformation(oriented_pose, rotation_matrix, position)
                
                orientations.append(oriented_pose)
                
            except Exception as e:
                self.logger.debug(f"Error generating orientation {i}: {e}")
                continue
        
        return orientations
    
    def _simulated_annealing(self, protein: Any, initial_pose: Any) -> Any:
        """Perform simulated annealing optimization."""
        current_pose = deepcopy(initial_pose)
        current_score = self.scoring_function.score(protein, current_pose)
        
        best_pose = deepcopy(current_pose)
        best_score = current_score
        
        # Initialize temperature
        temperature = self.high_temp
        
        accepted_moves = 0
        total_moves = 0
        
        for step in range(self.md_steps):
            try:
                # Generate new pose by small perturbation
                new_pose = deepcopy(current_pose)
                self._apply_sa_perturbation(new_pose, temperature)
                
                # Evaluate new pose
                new_score = self.scoring_function.score(protein, new_pose)
                
                # Metropolis acceptance criterion
                delta_e = new_score - current_score
                
                if delta_e < 0 or np.random.random() < np.exp(-delta_e / (self.kb * temperature)):
                    # Accept move
                    current_pose = new_pose
                    current_score = new_score
                    accepted_moves += 1
                    
                    # Update best if better
                    if new_score < best_score:
                        best_pose = deepcopy(new_pose)
                        best_score = new_score
                
                total_moves += 1
                
                # Cool down temperature
                if step % 50 == 0:
                    temperature = max(self.target_temp, temperature * self.cooling_factor)
                
            except Exception as e:
                self.logger.debug(f"Error in SA step {step}: {e}")
                continue
        
        acceptance_rate = accepted_moves / max(total_moves, 1)
        self.logger.debug(f"SA completed: {acceptance_rate:.2%} acceptance rate")
        
        return best_pose
    
    def _apply_sa_perturbation(self, pose: Any, temperature: float) -> None:
        """Apply small perturbation for simulated annealing."""
        try:
            # Temperature-dependent perturbation magnitude
            max_translation = 1.0 * (temperature / self.high_temp)  # Up to 1 Å
            max_rotation = 0.3 * (temperature / self.high_temp)     # Up to ~17 degrees
            
            # Random translation
            translation = np.random.normal(0, max_translation, 3)
            
            # Random rotation around centroid
            rotation_angles = np.random.normal(0, max_rotation, 3)
            rotation_matrix = rotation_matrix_from_euler(rotation_angles)
            
            # Get current centroid
            centroid = self._get_pose_centroid(pose)
            if centroid is not None:
                # Rotate around centroid
                self._rotate_around_point(pose, rotation_matrix, centroid)
            
            # Apply translation
            self._translate_pose(pose, translation)
            
        except Exception as e:
            self.logger.debug(f"Error applying SA perturbation: {e}")
    
    def _final_minimization(self, protein: Any, pose: Any) -> Any:
        """Perform final energy minimization."""
        try:
            best_pose = deepcopy(pose)
            best_score = self.scoring_function.score(protein, best_pose)
            
            # Gradient-free minimization using small random steps
            step_size = 0.1  # Start with small steps
            
            for step in range(self.minimize_steps):
                # Generate small random perturbation
                test_pose = deepcopy(best_pose)
                
                # Small random translation
                translation = np.random.normal(0, step_size, 3)
                self._translate_pose(test_pose, translation)
                
                # Small random rotation
                rotation_angles = np.random.normal(0, step_size * 0.1, 3)
                rotation_matrix = rotation_matrix_from_euler(rotation_angles)
                
                centroid = self._get_pose_centroid(test_pose)
                if centroid is not None:
                    self._rotate_around_point(test_pose, rotation_matrix, centroid)
                
                # Evaluate new pose
                new_score = self.scoring_function.score(protein, test_pose)
                
                # Accept if better
                if new_score < best_score:
                    best_pose = test_pose
                    best_score = new_score
                    step_size *= 1.1  # Increase step size on success
                else:
                    step_size *= 0.9  # Decrease step size on failure
                
                # Keep step size reasonable
                step_size = max(0.01, min(0.5, step_size))
            
            return best_pose
            
        except Exception as e:
            self.logger.debug(f"Error in final minimization: {e}")
            return pose
    
    def _check_pose_validity(self, protein: Any, pose: Any, score: float) -> bool:
        """Check if pose is valid using soft energy criteria."""
        try:
            # Soft energy threshold (less strict than clash detection)
            energy_threshold = 100.0  # kcal/mol
            
            if score > energy_threshold:
                return False
            
            # Basic clash check (less strict)
            return self._basic_clash_check(protein, pose)
            
        except Exception as e:
            self.logger.debug(f"Error in pose validity check: {e}")
            return True  # Default to valid
    
    def _basic_clash_check(self, protein: Any, pose: Any, 
                          clash_threshold: float = 1.5) -> bool:
        """Basic clash detection with soft threshold."""
        try:
            protein_atoms = self._get_protein_atoms(protein)
            pose_atoms = self._get_pose_atoms(pose)
            
            if not protein_atoms or not pose_atoms:
                return True
            
            # Count severe clashes
            clash_count = 0
            max_clashes = 5  # Allow some flexibility
            
            for p_atom in protein_atoms[:100]:  # Limit for performance
                p_coords = self._get_atom_coords(p_atom)
                if p_coords is None:
                    continue
                
                for l_atom in pose_atoms:
                    l_coords = self._get_atom_coords(l_atom)
                    if l_coords is None:
                        continue
                    
                    distance = np.linalg.norm(p_coords - l_coords)
                    
                    if distance < clash_threshold:
                        clash_count += 1
                        if clash_count > max_clashes:
                            return False
            
            return True
            
        except Exception as e:
            self.logger.debug(f"Error in clash check: {e}")
            return True
    
    # Helper methods (shared with RandomSearch)
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
                # Extract from atoms list
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
                # Set in atoms list
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
        """Get protein atoms for analysis."""
        try:
            if hasattr(protein, 'active_site') and protein.active_site and 'atoms' in protein.active_site:
                return protein.active_site['atoms']
            elif hasattr(protein, 'atoms'):
                return protein.atoms
            else:
                return []
        except:
            return []
    
    def _get_pose_atoms(self, pose: Any) -> List[Any]:
        """Get pose atoms for analysis."""
        try:
            if hasattr(pose, 'atoms'):
                return pose.atoms
            elif hasattr(pose, 'molecule') and hasattr(pose.molecule, 'atoms'):
                return pose.molecule.atoms
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