"""
Enhanced genetic algorithm implementation for molecular docking.

This module provides a comprehensive genetic algorithm with advanced features
including tournament selection, blend crossover, adaptive mutation, and
local optimization capabilities based on the original PandaDock implementation.
"""

import numpy as np
import logging
from typing import List, Tuple, Any, Dict, Optional
import random
from copy import deepcopy
import time

from .base_algorithm import BaseAlgorithm
from ..utils_new.math_utils import rotation_matrix_from_euler, apply_rotation_translation

logger = logging.getLogger(__name__)


class GeneticAlgorithm(BaseAlgorithm):
    """
    Advanced genetic algorithm for molecular docking.
    
    Features:
    - Tournament selection with configurable size
    - Blend crossover for smooth pose interpolation
    - Adaptive mutation with both small and large perturbations
    - Elitist selection maintaining best individuals
    - Optional local optimization of top poses
    - Convergence detection and early stopping
    - Search radius adaptation for focused exploration
    """
    
    def __init__(self, scoring_function: Any, population_size: int = 50,
                 max_generations: int = 100, mutation_rate: float = 0.1, 
                 crossover_rate: float = 0.8, tournament_size: int = 3,
                 elitism_count: int = 2, adaptive_mutation: bool = True,
                 local_optimization: bool = True, convergence_patience: int = 20,
                 initial_radius: float = 10.0, **kwargs):
        """
        Initialize genetic algorithm.
        
        Args:
            scoring_function: Function to evaluate poses
            population_size: Size of the population
            max_generations: Maximum number of generations
            mutation_rate: Initial probability of mutation
            crossover_rate: Probability of crossover
            tournament_size: Size of tournament for selection
            elitism_count: Number of best individuals to preserve
            adaptive_mutation: Whether to use adaptive mutation rates
            local_optimization: Whether to optimize top poses
            convergence_patience: Generations without improvement before stopping
            initial_radius: Initial search radius around active site
        """
        # Set max_iterations for base class
        kwargs['max_iterations'] = max_generations
        super().__init__(scoring_function, **kwargs)
        
        self.population_size = population_size
        self.max_generations = max_generations
        self.mutation_rate = mutation_rate
        self.crossover_rate = crossover_rate
        self.tournament_size = tournament_size
        self.elitism_count = elitism_count
        self.adaptive_mutation = adaptive_mutation
        self.local_optimization = local_optimization
        self.convergence_patience = convergence_patience
        self.initial_radius = initial_radius
        
        # Algorithm state
        self.population = []
        self.fitness_scores = []
        self.best_individual = None
        self.best_score = float('inf')
        self.generations_without_improvement = 0
        self.diversity_history = []
        
        # Adaptive parameters
        self.current_mutation_rate = mutation_rate
        self.search_radius = initial_radius
        
        self.logger = logging.getLogger(__name__)
        self.logger.info(
            f"Initialized genetic algorithm: pop_size={population_size}, "
            f"max_gen={max_generations}, mutation_rate={mutation_rate}, "
            f"tournament_size={tournament_size}, elitism={elitism_count}"
        )
    
    def search(self, protein: Any, ligand: Any) -> List[Tuple[Any, float]]:
        """
        Run genetic algorithm search.
        
        Args:
            protein: Protein molecule object
            ligand: Ligand molecule object
            
        Returns:
            List of (pose, score) tuples sorted by score
        """
        if not self.validate_inputs(protein, ligand):
            raise ValueError("Invalid protein or ligand inputs")
        
        self.logger.info("Starting genetic algorithm search")
        start_time = time.time()
        
        # Get active site information
        active_site = self._get_active_site(protein)
        if not active_site:
            raise ValueError("No active site defined for protein")
        
        # Initialize population
        self._initialize_population(protein, ligand, active_site)
        
        # Evolution loop
        for generation in range(self.max_generations):
            self.iteration_count = generation
            
            # Evaluate fitness
            self._evaluate_population(protein, ligand)
            
            # Calculate population diversity
            diversity = self._calculate_diversity()
            self.diversity_history.append(diversity)
            
            # Log progress
            if generation % 10 == 0 or generation < 5:
                best_idx = np.argmin(self.fitness_scores)
                avg_score = np.mean([s for s in self.fitness_scores if s != float('inf')])
                self.logger.info(
                    f"Gen {generation}: best={self.fitness_scores[best_idx]:.4f}, "
                    f"avg={avg_score:.4f}, diversity={diversity:.4f}, "
                    f"radius={self.search_radius:.2f}, mut_rate={self.current_mutation_rate:.3f}"
                )
            
            # Check convergence
            if self._check_convergence():
                self.logger.info(f"Converged at generation {generation}")
                break
            
            # Update adaptive parameters
            if self.adaptive_mutation:
                self._update_adaptive_parameters(generation)
            
            # Create next generation
            self._create_next_generation()
        
        elapsed = time.time() - start_time
        self.logger.info(f"Genetic algorithm completed in {elapsed:.1f}s")
        
        # Final evaluation
        self._evaluate_population(protein, ligand)
        
        # Apply local optimization if requested
        if self.local_optimization:
            self._apply_local_optimization(protein, ligand)
        
        # Prepare results
        results = []
        for i, (individual, score) in enumerate(zip(self.population, self.fitness_scores)):
            if score != float('inf'):
                pose = self._individual_to_pose(individual, ligand)
                results.append((pose, score))
        
        # Sort results by score
        results.sort(key=lambda x: x[1])
        
        self.logger.info(
            f"Genetic algorithm final results: {len(results)} valid poses, "
            f"best_score={results[0][1] if results else 'N/A':.4f}"
        )
        
        return results
    
    def _initialize_population(self, protein: Any, ligand: Any, active_site: Dict) -> None:
        """Initialize random population of poses within active site."""
        self.population = []
        self.fitness_scores = []
        
        center = active_site['center']
        radius = active_site.get('radius', self.initial_radius)
        self.search_radius = min(radius, self.initial_radius)
        
        # Get reference ligand coordinates
        ref_coords = self._get_ligand_coordinates(ligand)
        if ref_coords is None:
            raise ValueError("Could not extract ligand coordinates")
        
        for i in range(self.population_size):
            try:
                individual = self._generate_random_individual(ref_coords, center, self.search_radius)
                self.population.append(individual)
            except Exception as e:
                self.logger.warning(f"Failed to generate individual {i}: {e}")
                # Add a default individual as fallback
                self.population.append({
                    'translation': center + np.random.normal(0, 2, 3),
                    'rotation': np.random.uniform(0, 2*np.pi, 3),
                    'coords': ref_coords.copy()
                })
        
        self.logger.info(f"Initialized population of {len(self.population)} individuals")
    
    def _generate_random_individual(self, ref_coords: np.ndarray, 
                                  center: np.ndarray, radius: float) -> Dict:
        """Generate a random individual (pose) within the search space."""
        # Random position within sphere
        while True:
            random_point = np.random.uniform(-1, 1, 3)
            if np.linalg.norm(random_point) <= 1.0:
                break
        
        translation = center + random_point * radius * np.random.random()
        
        # Random orientation
        rotation = np.random.uniform(0, 2*np.pi, 3)
        
        # Apply transformation to coordinates
        rotation_matrix = rotation_matrix_from_euler(rotation)
        transformed_coords = apply_rotation_translation(ref_coords, rotation_matrix, translation)
        
        return {
            'translation': translation,
            'rotation': rotation,
            'coords': transformed_coords
        }
    
    def _evaluate_population(self, protein: Any, ligand: Any) -> None:
        """Evaluate fitness of all individuals in population."""
        new_scores = []
        
        for i, individual in enumerate(self.population):
            try:
                pose = self._individual_to_pose(individual, ligand)
                score = self.scoring_function.score(protein, pose)
                new_scores.append(score)
                
                # Track best individual
                if score < self.best_score:
                    self.best_score = score
                    self.best_individual = deepcopy(individual)
                    self.generations_without_improvement = 0
                    
                    # Shrink search radius when better poses are found
                    self.search_radius = max(2.0, self.search_radius * 0.95)
                    
            except Exception as e:
                self.logger.debug(f"Scoring failed for individual {i}: {e}")
                new_scores.append(float('inf'))
        
        self.fitness_scores = new_scores
        
        # Track improvement
        if len(new_scores) > 0 and min(new_scores) >= self.best_score:
            self.generations_without_improvement += 1
    
    def _individual_to_pose(self, individual: Dict, reference_ligand: Any) -> Any:
        """Convert individual representation to pose object."""
        try:
            pose = deepcopy(reference_ligand)
            
            # Apply coordinates from individual
            coords = individual['coords']
            
            # Set coordinates using various possible attributes
            if hasattr(pose, 'coords'):
                pose.coords = coords
            elif hasattr(pose, 'coordinates'):
                pose.coordinates = coords
            elif hasattr(pose, 'xyz'):
                pose.xyz = coords
            elif hasattr(pose, 'atoms'):
                # Update atom coordinates
                for i, atom in enumerate(pose.atoms):
                    if i < len(coords):
                        if isinstance(atom, dict):
                            atom['coords'] = coords[i]
                        elif hasattr(atom, 'coords'):
                            atom.coords = coords[i]
            
            return pose
            
        except Exception as e:
            self.logger.debug(f"Error converting individual to pose: {e}")
            return deepcopy(reference_ligand)
    
    def _calculate_diversity(self) -> float:
        """Calculate population diversity based on coordinate spread."""
        try:
            if len(self.population) < 2:
                return 0.0
            
            # Calculate centroid distances
            centroids = []
            for individual in self.population:
                coords = individual['coords']
                centroid = np.mean(coords, axis=0)
                centroids.append(centroid)
            
            centroids = np.array(centroids)
            
            # Calculate average pairwise distance
            total_distance = 0.0
            count = 0
            
            for i in range(len(centroids)):
                for j in range(i+1, len(centroids)):
                    distance = np.linalg.norm(centroids[i] - centroids[j])
                    total_distance += distance
                    count += 1
            
            if count > 0:
                return total_distance / count
            else:
                return 0.0
                
        except Exception as e:
            self.logger.debug(f"Error calculating diversity: {e}")
            return 0.0
    
    def _check_convergence(self) -> bool:
        """Check if algorithm has converged."""
        # Check patience
        if self.generations_without_improvement >= self.convergence_patience:
            return True
        
        # Check diversity
        if len(self.diversity_history) >= 10:
            recent_diversity = self.diversity_history[-10:]
            if np.mean(recent_diversity) < 1.0:  # Very low diversity
                return True
        
        return False
    
    def _update_adaptive_parameters(self, generation: int) -> None:
        """Update adaptive parameters based on search progress."""
        # Adaptive mutation rate
        if self.generations_without_improvement > 5:
            # Increase mutation rate if stuck
            self.current_mutation_rate = min(0.5, self.current_mutation_rate * 1.1)
        else:
            # Decrease mutation rate if improving
            self.current_mutation_rate = max(0.01, self.current_mutation_rate * 0.95)
        
        # Adaptive search radius
        if generation > 0 and generation % 20 == 0:
            # Expand radius occasionally to avoid local optima
            self.search_radius = min(self.initial_radius, self.search_radius * 1.2)
    
    def _create_next_generation(self) -> None:
        """Create next generation through selection, crossover, and mutation."""
        new_population = []
        
        # Elitism: keep best individuals
        if self.elitism_count > 0 and self.population:
            # Sort by fitness
            elite_indices = np.argsort(self.fitness_scores)[:self.elitism_count]
            for idx in elite_indices:
                if self.fitness_scores[idx] != float('inf'):
                    new_population.append(deepcopy(self.population[idx]))
        
        # Generate rest of population
        while len(new_population) < self.population_size:
            # Selection
            parent1 = self._tournament_selection()
            parent2 = self._tournament_selection()
            
            # Crossover
            if random.random() < self.crossover_rate:
                child1, child2 = self._crossover(parent1, parent2)
            else:
                child1, child2 = deepcopy(parent1), deepcopy(parent2)
            
            # Mutation
            if random.random() < self.current_mutation_rate:
                child1 = self._mutate(child1)
            if random.random() < self.current_mutation_rate:
                child2 = self._mutate(child2)
            
            new_population.extend([child1, child2])
        
        # Trim to exact population size
        self.population = new_population[:self.population_size]
    
    def _tournament_selection(self) -> Dict:
        """Select individual using tournament selection."""
        tournament_indices = random.sample(
            range(len(self.population)), 
            min(self.tournament_size, len(self.population))
        )
        
        # Find best in tournament (lowest score)
        best_idx = min(tournament_indices, key=lambda i: self.fitness_scores[i])
        
        return deepcopy(self.population[best_idx])
    
    def _crossover(self, parent1: Dict, parent2: Dict) -> Tuple[Dict, Dict]:
        """Perform blend crossover between two parents."""
        try:
            # Blend crossover parameter
            alpha = random.uniform(0.2, 0.8)
            
            # Create children
            child1 = deepcopy(parent1)
            child2 = deepcopy(parent2)
            
            # Blend translations
            child1['translation'] = alpha * parent1['translation'] + (1 - alpha) * parent2['translation']
            child2['translation'] = (1 - alpha) * parent1['translation'] + alpha * parent2['translation']
            
            # Blend rotations (careful with angles)
            child1['rotation'] = alpha * parent1['rotation'] + (1 - alpha) * parent2['rotation']
            child2['rotation'] = (1 - alpha) * parent1['rotation'] + alpha * parent2['rotation']
            
            # Update coordinates based on new transformations
            self._update_individual_coords(child1)
            self._update_individual_coords(child2)
            
            return child1, child2
            
        except Exception as e:
            self.logger.debug(f"Error in crossover: {e}")
            return deepcopy(parent1), deepcopy(parent2)
    
    def _mutate(self, individual: Dict) -> Dict:
        """Apply mutation to individual."""
        try:
            mutated = deepcopy(individual)
            
            # Small perturbations
            translation_noise = np.random.normal(0, 0.5, 3)
            rotation_noise = np.random.normal(0, 0.1, 3)
            
            mutated['translation'] += translation_noise
            mutated['rotation'] += rotation_noise
            
            # Occasionally apply larger mutations
            if random.random() < 0.1:
                # Large translation
                mutated['translation'] += np.random.normal(0, 2.0, 3)
            
            if random.random() < 0.1:
                # Large rotation
                mutated['rotation'] += np.random.uniform(-np.pi, np.pi, 3)
            
            # Update coordinates
            self._update_individual_coords(mutated)
            
            return mutated
            
        except Exception as e:
            self.logger.debug(f"Error in mutation: {e}")
            return deepcopy(individual)
    
    def _update_individual_coords(self, individual: Dict) -> None:
        """Update individual coordinates based on transformation parameters."""
        try:
            # Get reference coordinates (should be stored initially)
            if 'reference_coords' in individual:
                ref_coords = individual['reference_coords']
            else:
                # Use current coords as reference
                ref_coords = individual['coords']
            
            # Apply transformation
            rotation_matrix = rotation_matrix_from_euler(individual['rotation'])
            transformed_coords = apply_rotation_translation(
                ref_coords, rotation_matrix, individual['translation']
            )
            
            individual['coords'] = transformed_coords
            
        except Exception as e:
            self.logger.debug(f"Error updating individual coordinates: {e}")
    
    def _apply_local_optimization(self, protein: Any, ligand: Any) -> None:
        """Apply local optimization to top individuals."""
        try:
            self.logger.info("Applying local optimization to top poses")
            
            # Sort population by fitness
            sorted_indices = np.argsort(self.fitness_scores)
            n_optimize = min(5, len(self.population))  # Optimize top 5
            
            for i in range(n_optimize):
                idx = sorted_indices[i]
                if self.fitness_scores[idx] == float('inf'):
                    continue
                
                original_score = self.fitness_scores[idx]
                optimized_individual = self._local_minimize(
                    protein, ligand, self.population[idx], original_score
                )
                
                # Update if improved
                if optimized_individual is not None:
                    self.population[idx] = optimized_individual
            
            # Re-evaluate after optimization
            self._evaluate_population(protein, ligand)
            
        except Exception as e:
            self.logger.warning(f"Error in local optimization: {e}")
    
    def _local_minimize(self, protein: Any, ligand: Any, individual: Dict, 
                       initial_score: float) -> Optional[Dict]:
        """Perform local minimization of an individual."""
        try:
            best_individual = deepcopy(individual)
            best_score = initial_score
            
            step_size = 0.2  # Small step size for local optimization
            n_steps = 20
            
            for step in range(n_steps):
                # Small random perturbation
                test_individual = deepcopy(best_individual)
                
                # Perturb translation and rotation
                test_individual['translation'] += np.random.normal(0, step_size, 3)
                test_individual['rotation'] += np.random.normal(0, step_size * 0.1, 3)
                
                # Update coordinates
                self._update_individual_coords(test_individual)
                
                # Evaluate
                pose = self._individual_to_pose(test_individual, ligand)
                score = self.scoring_function.score(protein, pose)
                
                # Accept if better
                if score < best_score:
                    best_individual = test_individual
                    best_score = score
                    step_size *= 1.1  # Increase step size on success
                else:
                    step_size *= 0.9  # Decrease step size on failure
                
                # Keep step size reasonable
                step_size = max(0.05, min(1.0, step_size))
            
            return best_individual if best_score < initial_score else None
            
        except Exception as e:
            self.logger.debug(f"Error in local minimization: {e}")
            return None
    
    def _get_ligand_coordinates(self, ligand: Any) -> Optional[np.ndarray]:
        """Extract coordinates from ligand object."""
        try:
            if hasattr(ligand, 'coords'):
                coords = np.array(ligand.coords)
            elif hasattr(ligand, 'coordinates'):
                coords = np.array(ligand.coordinates)
            elif hasattr(ligand, 'xyz'):
                coords = np.array(ligand.xyz)
            elif hasattr(ligand, 'atoms'):
                coords_list = []
                for atom in ligand.atoms:
                    if isinstance(atom, dict) and 'coords' in atom:
                        coords_list.append(atom['coords'])
                    elif hasattr(atom, 'coords'):
                        coords_list.append(atom.coords)
                if coords_list:
                    coords = np.array(coords_list)
                else:
                    return None
            else:
                return None
            
            # Ensure 2D array
            if coords.ndim == 1:
                coords = coords.reshape(-1, 3)
            
            return coords
            
        except Exception as e:
            self.logger.debug(f"Error extracting ligand coordinates: {e}")
            return None
    
    def _get_active_site(self, protein: Any) -> Optional[Dict]:
        """Extract active site information from protein."""
        try:
            if hasattr(protein, 'active_site') and protein.active_site:
                return protein.active_site
            elif hasattr(protein, 'binding_site'):
                return protein.binding_site
            else:
                # Create default active site from protein center
                if hasattr(protein, 'atoms'):
                    coords = []
                    for atom in protein.atoms:
                        if isinstance(atom, dict) and 'coords' in atom:
                            coords.append(atom['coords'])
                        elif hasattr(atom, 'coords'):
                            coords.append(atom.coords)
                    
                    if coords:
                        coords_array = np.array(coords)
                        center = np.mean(coords_array, axis=0)
                        max_dist = np.max(np.linalg.norm(coords_array - center, axis=1))
                        
                        return {
                            'center': center,
                            'radius': max_dist + 5.0
                        }
                
                # Final fallback
                return {
                    'center': np.array([0.0, 0.0, 0.0]),
                    'radius': 15.0
                }
                
        except Exception as e:
            self.logger.debug(f"Error getting active site: {e}")
            return {
                'center': np.array([0.0, 0.0, 0.0]),
                'radius': 15.0
            }


# Keep the simple version for compatibility
SimpleGeneticAlgorithm = GeneticAlgorithm