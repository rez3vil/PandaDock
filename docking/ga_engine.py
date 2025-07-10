"""
Genetic Algorithm-based docking engine (AutoDock Vina-style)
Implements fast virtual screening using evolutionary algorithms
"""

import numpy as np
from typing import List, Dict, Any, Optional, Tuple
import logging
from concurrent.futures import ThreadPoolExecutor
import random

from .base_engine import DockingEngine, Pose
from ..scoring.scoring_functions import ScoringFunctions
from ..utils.math_utils import rotation_matrix, quaternion_to_matrix


class Individual:
    """Represents an individual in the genetic algorithm"""
    
    def __init__(self, genes: np.ndarray, fitness: float = float('inf')):
        self.genes = genes.copy()  # [x, y, z, qw, qx, qy, qz, torsion1, torsion2, ...]
        self.fitness = fitness
        self.pose = None
        self.evaluated = False
    
    def copy(self):
        """Create a copy of this individual"""
        new_individual = Individual(self.genes, self.fitness)
        new_individual.evaluated = self.evaluated
        if self.pose:
            new_individual.pose = self.pose
        return new_individual


class GAEngine(DockingEngine):
    """
    Genetic Algorithm-based docking engine similar to AutoDock Vina
    
    Features:
    - Fast evolutionary search
    - Parallel evaluation
    - Adaptive parameters
    - Empirical scoring
    - Virtual screening optimization
    """
    
    def __init__(self, config):
        super().__init__(config)
        self.scoring = ScoringFunctions(config)
        
        # GA parameters
        self.population_size = config.docking.population_size
        self.generations = config.docking.generations
        self.mutation_rate = config.docking.mutation_rate
        self.crossover_rate = config.docking.crossover_rate
        self.elitism_rate = config.docking.elitism_rate
        
        # Gene encoding
        self.num_position_genes = 3  # x, y, z
        self.num_orientation_genes = 4  # quaternion qw, qx, qy, qz
        self.num_torsion_genes = 0  # Will be set based on ligand
        self.gene_bounds = {}
        
        # Search parameters
        self.local_search_rate = 0.1  # Fraction of population for local search
        self.diversity_threshold = 0.8
        self.stagnation_limit = 100
        
        # Performance optimization
        self.use_parallel = config.n_jobs > 1
        self.n_jobs = config.n_jobs
        
        self.logger.info("Initialized GAEngine (Vina-style)")
    
    def dock(self, protein_file: str, ligand_file: str) -> List[Pose]:
        """
        Main docking method using genetic algorithm
        
        Steps:
        1. Prepare receptor and ligand
        2. Initialize population
        3. Evolve population over generations
        4. Apply local search to best individuals
        5. Return top poses
        """
        self.logger.info(f"Starting GA-based docking: {protein_file} + {ligand_file}")
        
        # Prepare structures
        self.prepare_receptor(protein_file)
        self.prepare_ligand(ligand_file)
        
        # Initialize gene encoding
        self._setup_gene_encoding()
        
        # Initialize population
        population = self._initialize_population()
        self.logger.info(f"Initialized population of {len(population)} individuals")
        
        # Evolve population
        best_individuals = self._evolve_population(population)
        
        # Convert best individuals to poses
        final_poses = []
        for individual in best_individuals:
            pose = self._individual_to_pose(individual)
            final_poses.append(pose)
        
        # Sort by score and take top poses
        final_poses.sort(key=lambda x: x.score)
        final_poses = final_poses[:self.config.docking.num_poses]
        
        self.logger.info(f"Final result: {len(final_poses)} poses")
        return final_poses
    
    def _setup_gene_encoding(self):
        """Setup gene encoding based on ligand structure"""
        # In real implementation, this would analyze ligand structure
        # and determine rotatable bonds
        self.num_torsion_genes = 6  # Placeholder
        
        # Set gene bounds
        min_bounds, max_bounds = self.grid_box.get_bounds()
        
        self.gene_bounds = {
            'position': (min_bounds, max_bounds),
            'orientation': (np.array([-1, -1, -1, -1]), np.array([1, 1, 1, 1])),
            'torsions': (np.array([-np.pi] * self.num_torsion_genes), 
                        np.array([np.pi] * self.num_torsion_genes))
        }
        
        self.total_genes = (self.num_position_genes + 
                           self.num_orientation_genes + 
                           self.num_torsion_genes)
    
    def _initialize_population(self) -> List[Individual]:
        """Initialize random population"""
        population = []
        
        for _ in range(self.population_size):
            genes = self._generate_random_genes()
            individual = Individual(genes)
            population.append(individual)
        
        return population
    
    def _generate_random_genes(self) -> np.ndarray:
        """Generate random genes within bounds"""
        genes = np.zeros(self.total_genes)
        
        # Position genes
        pos_min, pos_max = self.gene_bounds['position']
        genes[:3] = np.random.uniform(pos_min, pos_max)
        
        # Orientation genes (quaternion)
        quat = np.random.randn(4)
        quat = quat / np.linalg.norm(quat)  # Normalize
        genes[3:7] = quat
        
        # Torsion genes
        torsion_min, torsion_max = self.gene_bounds['torsions']
        genes[7:] = np.random.uniform(torsion_min, torsion_max)
        
        return genes
    
    def _evolve_population(self, population: List[Individual]) -> List[Individual]:
        """Evolve population using genetic algorithm"""
        best_fitness_history = []
        stagnation_count = 0
        
        for generation in range(self.generations):
            # Evaluate population
            self._evaluate_population(population)
            
            # Track best fitness
            best_fitness = min(individual.fitness for individual in population)
            best_fitness_history.append(best_fitness)
            
            # Check for stagnation
            if len(best_fitness_history) > self.stagnation_limit:
                recent_best = min(best_fitness_history[-self.stagnation_limit:])
                if abs(recent_best - best_fitness) < 0.001:
                    stagnation_count += 1
                else:
                    stagnation_count = 0
            
            # Log progress
            if generation % 1000 == 0:
                avg_fitness = np.mean([ind.fitness for ind in population])
                self.logger.info(f"Generation {generation}: Best={best_fitness:.3f}, Avg={avg_fitness:.3f}")
            
            # Early stopping if stagnated
            if stagnation_count > 50:
                self.logger.info(f"Early stopping at generation {generation} due to stagnation")
                break
            
            # Selection
            parents = self._tournament_selection(population)
            
            # Crossover and mutation
            offspring = self._create_offspring(parents)
            
            # Apply local search to some individuals
            if random.random() < self.local_search_rate:
                self._apply_local_search(offspring)
            
            # Replacement
            population = self._replacement(population, offspring)
        
        # Return best individuals
        population.sort(key=lambda x: x.fitness)
        return population[:self.config.docking.num_poses * 2]  # Return more for clustering
    
    def _evaluate_population(self, population: List[Individual]):
        """Evaluate fitness of all individuals in population"""
        unevaluated = [ind for ind in population if not ind.evaluated]
        
        if not unevaluated:
            return
        
        if self.use_parallel:
            self._evaluate_parallel(unevaluated)
        else:
            self._evaluate_sequential(unevaluated)
    
    def _evaluate_sequential(self, individuals: List[Individual]):
        """Evaluate individuals sequentially"""
        for individual in individuals:
            individual.fitness = self._evaluate_individual(individual)
            individual.evaluated = True
    
    def _evaluate_parallel(self, individuals: List[Individual]):
        """Evaluate individuals in parallel"""
        with ThreadPoolExecutor(max_workers=self.n_jobs) as executor:
            fitness_values = list(executor.map(self._evaluate_individual, individuals))
        
        for individual, fitness in zip(individuals, fitness_values):
            individual.fitness = fitness
            individual.evaluated = True
    
    def _evaluate_individual(self, individual: Individual) -> float:
        """Evaluate fitness of a single individual"""
        try:
            # Convert genes to pose
            pose = self._genes_to_pose(individual.genes)
            
            # Check if pose is valid
            if not self.validate_pose(pose):
                return 1000.0  # High penalty for invalid poses
            
            # Calculate fitness (lower is better)
            fitness = self.score(pose)
            
            # Store pose in individual
            individual.pose = pose
            
            return fitness
            
        except Exception as e:
            self.logger.warning(f"Error evaluating individual: {e}")
            return 1000.0
    
    def _genes_to_pose(self, genes: np.ndarray) -> Pose:
        """Convert genes to pose coordinates"""
        # Extract gene components
        position = genes[:3]
        quaternion = genes[3:7]
        torsions = genes[7:]
        
        # Normalize quaternion
        quaternion = quaternion / np.linalg.norm(quaternion)
        
        # Generate ligand coordinates from torsions
        coords = self._build_ligand_from_torsions(torsions)
        
        # Apply rotation
        rotation_matrix = quaternion_to_matrix(quaternion)
        rotated_coords = np.dot(coords, rotation_matrix.T)
        
        # Apply translation
        final_coords = rotated_coords + position
        
        # Create pose
        pose = Pose(
            coordinates=final_coords,
            score=0.0,
            energy=0.0,
            ligand_name=self.ligand,
            pose_id=f"ga_pose_{id(genes)}"
        )
        
        return pose
    
    def _build_ligand_from_torsions(self, torsions: np.ndarray) -> np.ndarray:
        """Build ligand coordinates from torsion angles"""
        # Placeholder implementation
        # In real implementation, this would use molecular mechanics
        # to build coordinates from torsion angles
        
        num_atoms = 20  # Placeholder
        coords = np.random.randn(num_atoms, 3) * 2.0
        
        # Apply torsion-dependent perturbations
        for i, torsion in enumerate(torsions):
            angle_factor = np.sin(torsion) * 0.5
            coords[i % num_atoms] += angle_factor
        
        return coords
    
    def _individual_to_pose(self, individual: Individual) -> Pose:
        """Convert individual to pose with final scoring"""
        if individual.pose is None:
            individual.pose = self._genes_to_pose(individual.genes)
        
        # Update score
        individual.pose.score = individual.fitness
        
        return individual.pose
    
    def _tournament_selection(self, population: List[Individual], tournament_size: int = 3) -> List[Individual]:
        """Select parents using tournament selection"""
        parents = []
        
        for _ in range(len(population)):
            # Select random individuals for tournament
            tournament = random.sample(population, tournament_size)
            
            # Select best from tournament
            winner = min(tournament, key=lambda x: x.fitness)
            parents.append(winner.copy())
        
        return parents
    
    def _create_offspring(self, parents: List[Individual]) -> List[Individual]:
        """Create offspring through crossover and mutation"""
        offspring = []
        
        for i in range(0, len(parents), 2):
            parent1 = parents[i]
            parent2 = parents[i + 1] if i + 1 < len(parents) else parents[0]
            
            # Crossover
            if random.random() < self.crossover_rate:
                child1, child2 = self._crossover(parent1, parent2)
            else:
                child1, child2 = parent1.copy(), parent2.copy()
            
            # Mutation
            if random.random() < self.mutation_rate:
                self._mutate(child1)
            if random.random() < self.mutation_rate:
                self._mutate(child2)
            
            offspring.extend([child1, child2])
        
        return offspring
    
    def _crossover(self, parent1: Individual, parent2: Individual) -> Tuple[Individual, Individual]:
        """Perform crossover between two parents"""
        # Single-point crossover
        crossover_point = random.randint(1, len(parent1.genes) - 1)
        
        child1_genes = np.concatenate([
            parent1.genes[:crossover_point],
            parent2.genes[crossover_point:]
        ])
        
        child2_genes = np.concatenate([
            parent2.genes[:crossover_point],
            parent1.genes[crossover_point:]
        ])
        
        # Ensure quaternion normalization
        child1_genes[3:7] = child1_genes[3:7] / np.linalg.norm(child1_genes[3:7])
        child2_genes[3:7] = child2_genes[3:7] / np.linalg.norm(child2_genes[3:7])
        
        return Individual(child1_genes), Individual(child2_genes)
    
    def _mutate(self, individual: Individual):
        """Mutate an individual"""
        genes = individual.genes.copy()
        
        # Position mutation
        if random.random() < 0.3:
            pos_noise = np.random.normal(0, 0.5, 3)
            genes[:3] += pos_noise
            
            # Ensure within bounds
            pos_min, pos_max = self.gene_bounds['position']
            genes[:3] = np.clip(genes[:3], pos_min, pos_max)
        
        # Orientation mutation
        if random.random() < 0.3:
            quat_noise = np.random.normal(0, 0.1, 4)
            genes[3:7] += quat_noise
            genes[3:7] = genes[3:7] / np.linalg.norm(genes[3:7])
        
        # Torsion mutation
        if random.random() < 0.5:
            torsion_indices = np.random.choice(self.num_torsion_genes, 
                                             size=random.randint(1, self.num_torsion_genes), 
                                             replace=False)
            for idx in torsion_indices:
                genes[7 + idx] += np.random.normal(0, 0.3)
                genes[7 + idx] = np.clip(genes[7 + idx], -np.pi, np.pi)
        
        individual.genes = genes
        individual.evaluated = False
    
    def _apply_local_search(self, individuals: List[Individual]):
        """Apply local search to improve individuals"""
        # Select random individuals for local search
        num_local_search = max(1, int(len(individuals) * self.local_search_rate))
        selected = random.sample(individuals, num_local_search)
        
        for individual in selected:
            self._local_search(individual)
    
    def _local_search(self, individual: Individual):
        """Perform local search on an individual"""
        best_fitness = individual.fitness
        best_genes = individual.genes.copy()
        
        # Try small perturbations
        for _ in range(10):
            # Create perturbed copy
            perturbed_genes = individual.genes.copy()
            
            # Add small random perturbations
            perturbed_genes[:3] += np.random.normal(0, 0.1, 3)  # Position
            perturbed_genes[3:7] += np.random.normal(0, 0.05, 4)  # Orientation
            perturbed_genes[3:7] = perturbed_genes[3:7] / np.linalg.norm(perturbed_genes[3:7])
            perturbed_genes[7:] += np.random.normal(0, 0.1, self.num_torsion_genes)  # Torsions
            
            # Evaluate
            perturbed_individual = Individual(perturbed_genes)
            fitness = self._evaluate_individual(perturbed_individual)
            
            # Keep if better
            if fitness < best_fitness:
                best_fitness = fitness
                best_genes = perturbed_genes.copy()
        
        # Update individual
        individual.genes = best_genes
        individual.fitness = best_fitness
        individual.evaluated = True
    
    def _replacement(self, population: List[Individual], offspring: List[Individual]) -> List[Individual]:
        """Replace population with offspring using elitism"""
        # Combine population and offspring
        combined = population + offspring
        
        # Sort by fitness
        combined.sort(key=lambda x: x.fitness)
        
        # Keep best individuals (elitism)
        new_population = combined[:self.population_size]
        
        return new_population
    
    def score(self, pose: Pose) -> float:
        """Score a pose using Vina-like scoring function"""
        # Use empirical scoring similar to Vina
        score = self.scoring.calculate_vina_score(pose.coordinates)
        
        # Add clash penalty
        clash_score = self.scoring.calculate_clash_score(pose.coordinates)
        score += clash_score * 5.0
        
        # Update pose scoring details
        pose.vdw_energy = self.scoring.calculate_vdw_energy(pose.coordinates)
        pose.hbond_energy = self.scoring.calculate_hbond_energy(pose.coordinates)
        pose.hydrophobic_energy = self.scoring.calculate_hydrophobic_energy(pose.coordinates)
        pose.clash_score = clash_score
        
        return score
    
    def get_engine_info(self) -> Dict[str, Any]:
        """Get information about the GA engine"""
        info = super().get_engine_info()
        info.update({
            'engine_type': 'GAEngine',
            'description': 'AutoDock Vina-style genetic algorithm docking',
            'features': [
                'Evolutionary search',
                'Parallel evaluation',
                'Local search optimization',
                'Empirical scoring',
                'Fast virtual screening'
            ],
            'parameters': {
                'population_size': self.population_size,
                'generations': self.generations,
                'mutation_rate': self.mutation_rate,
                'crossover_rate': self.crossover_rate,
                'elitism_rate': self.elitism_rate,
                'local_search_rate': self.local_search_rate,
                'n_jobs': self.n_jobs
            },
            'gene_encoding': {
                'position_genes': self.num_position_genes,
                'orientation_genes': self.num_orientation_genes,
                'torsion_genes': self.num_torsion_genes,
                'total_genes': self.total_genes
            }
        })
        return info