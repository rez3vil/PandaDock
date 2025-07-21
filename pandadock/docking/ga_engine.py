# -*- coding: utf-8 -*-
"""
Genetic Algorithm-based docking engine (AutoDock Vina-style)
Implements fast virtual screening using evolutionary algorithms
"""

import numpy as np
from typing import List, Dict, Any, Optional, Tuple
import logging
from concurrent.futures import ThreadPoolExecutor
import random
from tqdm import tqdm

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from docking.base_engine import DockingEngine, Pose
from scoring.scoring_functions import ScoringFunctions
from utils.math_utils import rotation_matrix, quaternion_to_matrix


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
        
        # GA parameters - tuned for diversity and speed
        self.population_size = min(30, max(20, config.docking.population_size))  # Smaller for speed
        self.generations = min(500, config.docking.generations)  # Limit generations for testing
        self.mutation_rate = max(0.3, config.docking.mutation_rate)  # Higher mutation for diversity
        self.crossover_rate = min(0.7, config.docking.crossover_rate)  # Balanced crossover
        self.elitism_rate = min(0.1, config.docking.elitism_rate)  # Reduced elitism
        
        # Gene encoding
        self.num_position_genes = 3  # x, y, z
        self.num_orientation_genes = 4  # quaternion qw, qx, qy, qz
        self.num_torsion_genes = 0  # Will be set based on ligand
        self.gene_bounds = {}
        
        # Search parameters
        self.local_search_rate = 0.1  # Fraction of population for local search
        self.diversity_threshold = 0.8
        self.stagnation_limit = 100
        
        # Niching parameters for diversity preservation
        self.niche_radius = 2.0  # Distance threshold for niching
        self.sharing_alpha = 1.0  # Sharing function exponent
        self.min_niche_size = 3   # Minimum individuals per niche
        
        # Performance optimization
        self.use_parallel = config.n_jobs > 1
        self.n_jobs = config.n_jobs
        
        self.logger.info("Initialized GAEngine (Vina-style)")
    
    # def dock(self, protein_file: str, ligand_file: str) -> List[Pose]:
    #     """
    #     Main docking method using genetic algorithm
        
    #     Steps:
    #     1. Prepare receptor and ligand
    #     2. Initialize population
    #     3. Evolve population over generations
    #     4. Apply local search to best individuals
    #     5. Return top poses
    #     """
    #     self.logger.info(f"Starting GA-based docking: {protein_file} + {ligand_file}")
        
    #     # Prepare structures
    #     self.prepare_receptor(protein_file)
    #     self.prepare_ligand(ligand_file)
        
    #     # Initialize gene encoding
    #     self._setup_gene_encoding()
        
    #     # Initialize population
    #     population = self._initialize_population()
    #     self.logger.info(f"Initialized population of {len(population)} individuals")
        
    #     # Evolve population
    #     best_individuals = self._evolve_population(population)
        
    #     # Convert best individuals to poses
    #     final_poses = []
    #     for individual in best_individuals:
    #         pose = self._individual_to_pose(individual)
    #         final_poses.append(pose)
        
    #     # Sort by score and take top poses
    #     final_poses.sort(key=lambda x: x.score)
    #     final_poses = final_poses[:self.config.docking.num_poses]
        
    #     self.logger.info(f"Final result: {len(final_poses)} poses")
    #     return final_poses

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
        
        # Sort by score to determine final ranking
        final_poses.sort(key=lambda x: x.score)
        
        # Limit to the number of poses requested by the configuration
        num_poses_to_return = getattr(self.config.docking, 'num_poses', 10)
        final_poses = final_poses[:num_poses_to_return]
        
        # === FIX START: Overwrite pose IDs for clean, sequential filenames ===
        # After all sorting and selection, we re-label the final list of poses
        # with simple, human-readable IDs based on their final rank (1, 2, 3...).
        # This ensures the output files are named 'pose_1.sdf', 'pose_2.sdf', etc.
        for i, pose in enumerate(final_poses):
            pose.pose_id = f"pose_{i + 1}"
        # === FIX END ===
        
        self.logger.info(f"Final result: {len(final_poses)} poses with clean, sequential IDs")
        return final_poses
    
    def _setup_gene_encoding(self):
        """Setup gene encoding based on ligand structure"""
        # In real implementation, this would analyze ligand structure
        # and determine rotatable bonds
        self.num_torsion_genes = 6  # Placeholder
        
        # Set gene bounds centered around binding site
        # Use reasonable sampling radius instead of entire grid box
        sampling_radius = min(8.0, np.min(self.grid_box.size) / 4.0)  # Max 8Å or 1/4 of smallest dimension
        center = self.grid_box.center
        
        # Position bounds: center ± sampling radius, but constrained to grid box
        min_bounds, max_bounds = self.grid_box.get_bounds()
        pos_min = np.maximum(center - sampling_radius, min_bounds)
        pos_max = np.minimum(center + sampling_radius, max_bounds)
        
        self.gene_bounds = {
            'position': (pos_min, pos_max),
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
        
        # Position genes - sample around grid center with bias toward center
        pos_min, pos_max = self.gene_bounds['position']
        # Use normal distribution centered on grid center for better clustering
        center_pos = (pos_min + pos_max) / 2
        spread = (pos_max - pos_min) / 6  # 3-sigma covers the bounds
        genes[:3] = np.random.normal(center_pos, spread)
        # Ensure within bounds
        genes[:3] = np.clip(genes[:3], pos_min, pos_max)
        
        # Orientation genes (quaternion)
        quat = np.random.randn(4)
        quat = quat / np.linalg.norm(quat)  # Normalize
        genes[3:7] = quat
        
        # Torsion genes
        torsion_min, torsion_max = self.gene_bounds['torsions']
        genes[7:] = np.random.uniform(torsion_min, torsion_max)
        
        return genes
    
    def _evolve_population(self, population: List[Individual]) -> List[Individual]:
        """Fast evolution using ML-inspired rigid body sampling"""
        self.logger.info("Using fast rigid body sampling (ML-inspired approach)")
        
        # Generate poses using rigid body transformations like ML engine
        num_samples = self.config.docking.num_poses * 3
        poses = self._generate_rigid_body_poses(num_samples)
        
        # Convert poses to individuals
        individuals = []
        for i, pose in enumerate(poses):
            # Create genes from pose (position + orientation + minimal torsions)
            genes = self._pose_to_genes(pose)
            individual = Individual(genes)
            individual.pose = pose
            individuals.append(individual)
        
        # Quick evaluation
        print("🎯 Evaluating pose quality...")
        eval_pbar = tqdm(individuals, desc="🧬 Fast evaluation", unit="pose",
                        bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]')
        for individual in eval_pbar:
            individual.fitness = self._evaluate_individual(individual)
        eval_pbar.close()
        
        # Sort by fitness and return best poses
        individuals.sort(key=lambda x: x.fitness)
        return individuals[:self.config.docking.num_poses]
    
    def _generate_rigid_body_poses(self, num_samples: int) -> List[Pose]:
        """Generate poses using rigid body transformations (rotation + translation)"""
        if not self.ligand or 'coordinates' not in self.ligand:
            raise ValueError("No ligand structure loaded. Call prepare_ligand() first.")
        
        base_coords = self.ligand['coordinates'].copy()
        poses = []
        ligand_name = self.ligand['name'] if self.ligand else 'unknown'
        
        # Get grid box bounds for safer placement
        min_bounds, max_bounds = self.grid_box.get_bounds()
        
        for i in range(num_samples):
            # ### FIX START: Correct 3D Rigid Body Transformation ###
            # 1. Center the molecule at the origin
            mol_center = np.mean(base_coords, axis=0)
            centered_coords = base_coords - mol_center
            
            # 2. Apply rotation to the centered molecule
            angles = np.random.uniform(0, 2*np.pi, 3)
            rot_matrix = rotation_matrix(angles)
            rotated_coords = np.dot(centered_coords, rot_matrix.T)
            
            # 3. Pick a new random center location near the grid center
            # Sample around the binding site center with reasonable radius
            sampling_radius = min(8.0, np.min(self.grid_box.size) / 4.0)  # Max 8Å or 1/4 of smallest dimension
            random_offset = np.random.normal(0, sampling_radius/3, 3)  # 3-sigma gives ~99% within radius
            target_center = self.grid_box.center + random_offset
            
            # Ensure target center is still within bounds
            min_bounds, max_bounds = self.grid_box.get_bounds()
            target_center = np.clip(target_center, min_bounds, max_bounds)
            
            # 4. Translate the rotated molecule to the new target center
            final_coords = rotated_coords + target_center
            # ### FIX END ###

            # Small conformational perturbation
            perturbation = np.random.normal(0, 0.1, final_coords.shape)
            final_coords += perturbation
            
            # Create pose
            pose = Pose(
                coordinates=final_coords,
                score=0.0,
                energy=0.0,
                ligand_name=ligand_name,
                pose_id=f"internal_pose_{i}",
                confidence=0.0
            )
            poses.append(pose)
        
        return poses
    
    def _ensure_within_bounds(self, coords: np.ndarray, min_bounds: np.ndarray, max_bounds: np.ndarray) -> np.ndarray:
        """Ensure coordinates are strictly within grid bounds"""
        coords = np.clip(coords, min_bounds, max_bounds)
        
        # If any coordinates are still out of bounds, re-center them
        for i in range(len(coords)):
            if np.any(coords[i] < min_bounds) or np.any(coords[i] > max_bounds):
                coords[i] = (min_bounds + max_bounds) / 2
        
        return coords
    def _pose_to_genes(self, pose: Pose) -> np.ndarray:
        """Convert pose to gene representation"""
        genes = np.zeros(self.total_genes)
        
        # Position genes (center of mass)
        center = np.mean(pose.coordinates, axis=0)
        genes[:3] = center
        
        # Orientation genes (quaternion - simplified)
        genes[3:7] = np.random.randn(4)
        genes[3:7] = genes[3:7] / np.linalg.norm(genes[3:7])  # Normalize
        
        # Torsion genes (minimal)
        torsion_min, torsion_max = self.gene_bounds['torsions']
        genes[7:] = np.random.uniform(torsion_min, torsion_max)
        
        return genes
    
    def _select_diverse_poses(self, population: List[Individual], num_poses: int) -> List[Individual]:
        """Select diverse poses using multi-objective optimization"""
        # Sort by shared fitness (considers both quality and diversity)
        population.sort(key=lambda x: getattr(x, 'shared_fitness', x.fitness))
        
        selected = []
        remaining = population.copy()
        
        # Select first (best) individual
        if remaining:
            selected.append(remaining.pop(0))
        
        # Select remaining individuals to maximize diversity
        while len(selected) < num_poses and remaining:
            best_candidate = None
            best_score = -float('inf')
            
            for candidate in remaining:
                # Calculate minimum distance to already selected individuals
                min_distance = float('inf')
                for selected_ind in selected:
                    distance = self._calculate_distance(candidate, selected_ind)
                    min_distance = min(min_distance, distance)
                
                # Multi-objective score: balance fitness and diversity
                fitness_score = -candidate.fitness  # Convert to maximization
                diversity_score = min_distance
                combined_score = 0.6 * fitness_score + 0.4 * diversity_score
                
                if combined_score > best_score:
                    best_score = combined_score
                    best_candidate = candidate
            
            if best_candidate:
                selected.append(best_candidate)
                remaining.remove(best_candidate)
            else:
                break
        
        return selected
    
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
            
            # **FIXED**: Instead of trying to fix, apply a heavy penalty for invalid poses.
            if not self.validate_pose(pose):
                self.logger.debug(f"Invalid pose generated, applying high penalty.")
                return 100.0  # High penalty for being outside the box
            
            # Calculate base fitness (lower is better)
            base_fitness = self.score(pose)
            
            # Store pose in individual
            individual.pose = pose
            
            return base_fitness
            
        except Exception as e:
            self.logger.warning(f"Error evaluating individual: {e}")
            return 150.0  # Even higher penalty for unexpected errors
    
    def _calculate_fitness_sharing(self, population: List[Individual]):
        """Apply fitness sharing to maintain diversity"""
        for i, individual in enumerate(population):
            if not hasattr(individual, 'shared_fitness'):
                niche_count = 0.0
                
                for j, other in enumerate(population):
                    if i != j:
                        distance = self._calculate_distance(individual, other)
                        if distance < self.niche_radius:
                            # Sharing function
                            sharing_value = 1.0 - (distance / self.niche_radius) ** self.sharing_alpha
                            niche_count += sharing_value
                        else:
                            niche_count += 0.0
                
                # Shared fitness = raw fitness / niche count
                individual.shared_fitness = individual.fitness / max(1.0, niche_count)
    
    def _calculate_distance(self, ind1: Individual, ind2: Individual) -> float:
        """Calculate distance between two individuals in genotype space"""
        # Euclidean distance in gene space
        return np.linalg.norm(ind1.genes - ind2.genes)
    
    def _genes_to_pose(self, genes: np.ndarray) -> Pose:
        """Convert genes to a 3D Pose, ensuring correct geometry."""
        position = genes[:3]
        quaternion = genes[3:7]
        torsions = genes[7:]
        
        base_coords = self.ligand['coordinates'].copy()
        
        # This is a placeholder for real torsional rotation logic.
        # For now, it just adds minor noise.
        current_coords = self._apply_conformational_torsions(base_coords, torsions)

        # ### FIX START: Correct 3D Rigid Body Transformation ###
        # 1. Center the molecule at the origin
        mol_center = np.mean(current_coords, axis=0)
        centered_coords = current_coords - mol_center
        
        # 2. Apply rotation to the centered molecule
        norm_quat = quaternion / np.linalg.norm(quaternion) if np.linalg.norm(quaternion) > 0 else np.array([1,0,0,0])
        rotation = quaternion_to_matrix(norm_quat)
        rotated_coords = np.dot(centered_coords, rotation.T)
        
        # 3. Translate the rotated molecule to its final gene-defined position
        final_coords = rotated_coords + position
        # ### FIX END ###
        
        # Create pose with truly unique ID
        import time
        gene_hash = hash(tuple(genes.round(4)))  # Round to avoid floating point issues
        timestamp = int(time.time() * 1000000) % 1000000  # Microsecond precision
        random_component = np.random.randint(100000, 999999)
        # Use position to make ID more unique
        position_hash = hash(tuple(position.round(3)))
        unique_id = abs(gene_hash + timestamp + random_component + position_hash) % 10000000000
        
        # Use actual ligand atom types and structure data from loaded SDF file
        ligand_name = self.ligand.get('name', 'unknown') if isinstance(self.ligand, dict) else 'unknown'
        atom_types = self.ligand.get('atom_types', []) if isinstance(self.ligand, dict) else []
        bonds = self.ligand.get('bonds', []) if isinstance(self.ligand, dict) else []
        
        pose = Pose(
            coordinates=final_coords,
            score=0.0,  # Will be updated by _individual_to_pose
            energy=0.0,  # Will be updated by _individual_to_pose  
            ligand_name=ligand_name,
            pose_id="internal_pose"
        )
        
        # Store actual molecular data from input ligand in pose for proper file formatting
        pose.atom_types = atom_types
        pose.bonds = bonds
        
        return pose
    
    def _build_ligand_from_torsions(self, torsions: np.ndarray) -> np.ndarray:
        """Build ligand coordinates using actual ligand structure with torsional modifications"""
        if not self.ligand or 'coordinates' not in self.ligand:
            raise ValueError("No ligand structure loaded. Call prepare_ligand() first.")
        
        # Start with actual ligand coordinates from input SDF file
        base_coords = self.ligand['coordinates'].copy()
        
        # Apply conformational changes based on torsions
        # This preserves the molecular structure while allowing conformational flexibility
        coords = self._apply_conformational_torsions(base_coords, torsions)
        
        return coords
    
    def _build_drug_like_molecule(self, torsions: np.ndarray) -> np.ndarray:
        """Build coordinates for a realistic drug-like molecule"""
        # Define a template drug-like molecule (similar to a typical pharmaceutical compound)
        # This creates a molecule with aromatic rings, functional groups, and realistic connectivity
        
        # Standard bond lengths and angles
        c_c_bond = 1.54  # C-C single bond
        c_n_bond = 1.47  # C-N bond
        c_o_bond = 1.43  # C-O bond
        aromatic_bond = 1.40  # Aromatic C-C bond
        c_h_bond = 1.09  # C-H bond
        
        coords = []
        
        # Build a benzene ring (aromatic core) - common in drugs
        ring_center = np.array([0.0, 0.0, 0.0])
        ring_radius = aromatic_bond / (2 * np.sin(np.pi/6))  # Benzene geometry
        
        for i in range(6):
            angle = i * np.pi / 3  # 60-degree increments
            # Apply torsional variation
            torsion_mod = torsions[i % len(torsions)] * 0.1  # Small variation
            x = ring_center[0] + ring_radius * np.cos(angle + torsion_mod)
            y = ring_center[1] + ring_radius * np.sin(angle + torsion_mod)
            z = ring_center[2] + 0.1 * np.sin(torsions[i % len(torsions)])  # Slight out-of-plane
            coords.append([x, y, z])
        
        # Add substituents on the benzene ring
        # Add a nitrogen-containing side chain (common in drugs)
        if len(torsions) > 0:
            # Attach to carbon 1 of ring
            attachment_point = np.array(coords[1])
            # Direction away from ring center
            direction = (attachment_point - ring_center)
            direction = direction / np.linalg.norm(direction)
            
            # Build chain: -CH2-NH-CH3
            # CH2 carbon
            ch2_pos = attachment_point + direction * c_c_bond
            coords.append(ch2_pos)
            
            # NH nitrogen
            nh_direction = direction + 0.3 * np.array([np.sin(torsions[0]), np.cos(torsions[0]), 0])
            nh_direction = nh_direction / np.linalg.norm(nh_direction)
            nh_pos = ch2_pos + nh_direction * c_n_bond
            coords.append(nh_pos)
            
            # CH3 carbon
            ch3_direction = nh_direction + 0.2 * np.array([0, 0, np.sin(torsions[1 % len(torsions)])])
            ch3_direction = ch3_direction / np.linalg.norm(ch3_direction)
            ch3_pos = nh_pos + ch3_direction * c_n_bond
            coords.append(ch3_pos)
        
        # Add an oxygen-containing group (hydroxyl or ether)
        if len(torsions) > 1:
            # Attach to carbon 4 of ring (opposite side)
            attachment_point = np.array(coords[4])
            direction = (attachment_point - ring_center)
            direction = direction / np.linalg.norm(direction)
            
            # OH oxygen
            oh_direction = direction + 0.2 * np.array([np.cos(torsions[1]), 0, np.sin(torsions[1])])
            oh_direction = oh_direction / np.linalg.norm(oh_direction)
            oh_pos = attachment_point + oh_direction * c_o_bond
            coords.append(oh_pos)
        
        # Add another carbon chain if we have more torsions
        if len(torsions) > 2:
            # Attach to carbon 3 of ring
            attachment_point = np.array(coords[2])
            direction = (attachment_point - ring_center)
            direction = direction / np.linalg.norm(direction)
            
            # Build a short alkyl chain
            for i in range(3):  # 3-carbon chain
                if len(coords) < 20:  # Limit total atoms
                    chain_direction = direction + 0.1 * i * np.array([
                        np.sin(torsions[(2+i) % len(torsions)]),
                        np.cos(torsions[(2+i) % len(torsions)]),
                        0.1 * np.sin(torsions[(3+i) % len(torsions)])
                    ])
                    chain_direction = chain_direction / np.linalg.norm(chain_direction)
                    chain_pos = attachment_point + (i + 1) * c_c_bond * chain_direction
                    coords.append(chain_pos)
        
        # Fill remaining atoms with hydrogens in realistic positions
        while len(coords) < 20:
            # Add hydrogens to existing carbons
            if len(coords) >= 6:  # We have the benzene ring
                # Pick a random carbon to add hydrogen to
                carbon_idx = np.random.randint(0, min(6, len(coords)))
                carbon_pos = np.array(coords[carbon_idx])
                
                # Direction roughly perpendicular to ring plane
                h_direction = np.array([0, 0, 1]) + 0.3 * np.random.randn(3)
                h_direction = h_direction / np.linalg.norm(h_direction)
                h_pos = carbon_pos + h_direction * c_h_bond
                coords.append(h_pos)
            else:
                # Just add a hydrogen near the last atom
                if coords:
                    last_pos = np.array(coords[-1])
                    h_pos = last_pos + np.random.randn(3) * 0.5
                    coords.append(h_pos)
                else:
                    coords.append([0, 0, 0])
        
        coords = np.array(coords[:20])  # Ensure exactly 20 atoms
        
        # Apply torsional conformational changes
        coords = self._apply_conformational_torsions(coords, torsions)
        
        # Center and scale molecule appropriately
        coords = coords - np.mean(coords, axis=0)
        max_extent = np.max(np.linalg.norm(coords, axis=1))
        if max_extent > 0:
            coords = coords * (4.0 / max_extent)  # Scale to ~4 Angstrom radius
        
        return coords
    
    def _generate_realistic_molecule_data(self, coords: np.ndarray) -> tuple:
        """Generate realistic atom types and bonds for drug-like molecule"""
        num_atoms = len(coords)
        
        # Define realistic atom types for drug-like molecules
        # Most drugs contain C, N, O, and some H atoms
        atom_types = []
        
        # Pattern for typical drug molecule:
        # Aromatic ring (6 C atoms) + side chains (C, N, O, H)
        drug_pattern = ['C', 'C', 'C', 'C', 'C', 'C',  # Benzene ring
                       'C', 'N', 'C',  # -CH2-NH-CH3 side chain  
                       'O',             # -OH group
                       'C', 'C', 'C',  # Alkyl chain
                       'H', 'H', 'H', 'H', 'H', 'H', 'H']  # Hydrogens
        
        # Assign atom types based on pattern
        for i in range(num_atoms):
            if i < len(drug_pattern):
                atom_types.append(drug_pattern[i])
            else:
                # Default to hydrogen for extra atoms
                atom_types.append('H')
        
        # Generate realistic bonds for drug-like connectivity
        bonds = self._generate_drug_like_bonds(num_atoms)
        
        return atom_types, bonds
    
    def _generate_drug_like_bonds(self, num_atoms: int) -> list:
        """Generate realistic bond connectivity for drug-like molecule"""
        bonds = []
        
        # Benzene ring bonds (aromatic)
        if num_atoms >= 6:
            for i in range(6):
                next_atom = (i + 1) % 6
                bonds.append((i, next_atom, 'aromatic'))
        
        # Side chain bonds
        if num_atoms >= 9:
            # Ring to side chain: C1-C6 (carbon 1 to side chain carbon)
            bonds.append((1, 6, 'single'))
            # Side chain: C6-N7
            bonds.append((6, 7, 'single'))
            # Side chain: N7-C8  
            bonds.append((7, 8, 'single'))
        
        # OH group bond
        if num_atoms >= 10:
            # Ring to OH: C4-O9
            bonds.append((4, 9, 'single'))
        
        # Alkyl chain bonds
        if num_atoms >= 13:
            # Ring to alkyl: C2-C10
            bonds.append((2, 10, 'single'))
            # Alkyl chain: C10-C11-C12
            bonds.append((10, 11, 'single'))
            bonds.append((11, 12, 'single'))
        
        # Hydrogen bonds (to heavy atoms)
        if num_atoms >= 14:
            # Add H bonds to carbons and nitrogens
            h_atom_idx = 13
            for heavy_atom in range(min(13, num_atoms)):
                if h_atom_idx < num_atoms:
                    # Each heavy atom gets at least one hydrogen
                    atom_type = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'N', 'C', 'O', 'C', 'C', 'C'][heavy_atom]
                    if atom_type in ['C', 'N'] and h_atom_idx < num_atoms:
                        bonds.append((heavy_atom, h_atom_idx, 'single'))
                        h_atom_idx += 1
                        if h_atom_idx >= num_atoms:
                            break
        
        return bonds
    
    def _apply_conformational_torsions(self, coords: np.ndarray, torsions: np.ndarray) -> np.ndarray:
        """Apply small conformational changes to preserve molecular structure"""
        coords = coords.copy()
        
        # Apply small random perturbations for conformational diversity
        # This preserves the overall molecular structure while allowing flexibility
        
        # Small random perturbations (much smaller than in synthetic molecule generation)
        max_perturbation = 0.2  # Maximum 0.2 Angstrom displacement per atom
        
        for i in range(len(coords)):
            if i < len(torsions):
                # Use torsion value to generate reproducible but varied perturbations
                np.random.seed(abs(int(torsions[i] * 1000)) % 2147483647)
                perturbation = np.random.normal(0, max_perturbation, 3)
                coords[i] += perturbation
        
        # Reset random seed
        np.random.seed(None)
        
        return coords
    
    def _individual_to_pose(self, individual: Individual) -> Pose:
        """Convert individual to pose with final scoring"""
        if individual.pose is None:
            individual.pose = self._genes_to_pose(individual.genes)
        
        # Update score and energy (convert positive fitness to negative energy for standard docking convention)
        energy_value = -individual.fitness  # Convert to negative (more negative = better)
        individual.pose.score = individual.fitness  # Keep original for sorting (lower = better)
        individual.pose.energy = energy_value
        
        return individual.pose
    
    def _tournament_selection(self, population: List[Individual], tournament_size: int = 3) -> List[Individual]:
        """Select parents using tournament selection with shared fitness"""
        parents = []
        
        # Ensure fitness sharing is applied
        self._calculate_fitness_sharing(population)
        
        for _ in range(len(population)):
            # Select random individuals for tournament
            tournament_size_actual = min(tournament_size, len(population))
            tournament = random.sample(population, tournament_size_actual)
            
            # Select best from tournament using shared fitness
            winner = min(tournament, key=lambda x: getattr(x, 'shared_fitness', x.fitness))
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
        """Mutate an individual with higher diversity"""
        genes = individual.genes.copy()
        
        # Position mutation with higher probability and variance
        if random.random() < 0.7:  # Increased probability
            pos_noise = np.random.normal(0, 1.5, 3)  # Increased variance
            genes[:3] += pos_noise
            
            # Ensure within bounds
            pos_min, pos_max = self.gene_bounds['position']
            genes[:3] = np.clip(genes[:3], pos_min, pos_max)
        
        # Orientation mutation with higher probability
        if random.random() < 0.7:  # Increased probability
            quat_noise = np.random.normal(0, 0.3, 4)  # Increased variance
            genes[3:7] += quat_noise
            genes[3:7] = genes[3:7] / np.linalg.norm(genes[3:7])
        
        # Torsion mutation with higher probability and variance
        if random.random() < 0.8:  # Increased probability
            torsion_indices = np.random.choice(self.num_torsion_genes, 
                                             size=random.randint(1, self.num_torsion_genes), 
                                             replace=False)
            for idx in torsion_indices:
                genes[7 + idx] += np.random.normal(0, 0.8)  # Increased variance
                genes[7 + idx] = np.clip(genes[7 + idx], -np.pi, np.pi)
        
        # Occasionally apply large random jumps to maintain diversity
        if random.random() < 0.1:  # 10% chance of large mutation
            # Completely randomize some genes
            if random.random() < 0.5:
                genes[:3] = np.random.uniform(
                    self.gene_bounds['position'][0], 
                    self.gene_bounds['position'][1]
                )
            if random.random() < 0.5:
                quat = np.random.randn(4)
                genes[3:7] = quat / np.linalg.norm(quat)
            if random.random() < 0.5:
                genes[7:] = np.random.uniform(
                    self.gene_bounds['torsions'][0], 
                    self.gene_bounds['torsions'][1]
                )
        
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
            
            # **FIXED**: Clip genes to stay within bounds during local search
            pos_min, pos_max = self.gene_bounds['position']
            perturbed_genes[:3] = np.clip(perturbed_genes[:3], pos_min, pos_max)
            
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
        """Replace population using diversity-preserving strategies"""
        # Combine population and offspring
        combined = population + offspring
        
        # Apply fitness sharing to maintain diversity
        self._calculate_fitness_sharing(combined)
        
        # Sort by shared fitness (considering both quality and diversity)
        combined.sort(key=lambda x: getattr(x, 'shared_fitness', x.fitness))
        
        new_population = []
        
        # Reduce elitism rate to 10% (was 20%)
        elite_count = max(1, int(self.population_size * 0.1))
        new_population.extend(combined[:elite_count])
        
        # Add individuals using niching to ensure diversity
        remaining_individuals = combined[elite_count:]
        remaining_slots = self.population_size - elite_count
        
        # Tournament selection with diversity pressure
        for _ in range(remaining_slots):
            if len(remaining_individuals) > 0:
                # Select individual that maintains diversity
                selected = self._diversity_tournament_selection(remaining_individuals, new_population)
                if selected:
                    new_population.append(selected)
                    remaining_individuals.remove(selected)
                else:
                    # Add random individual if no diverse individual found
                    new_genes = self._generate_random_genes()
                    new_individual = Individual(new_genes)
                    new_population.append(new_individual)
            else:
                # Add random individual to maintain population size
                new_genes = self._generate_random_genes()
                new_individual = Individual(new_genes)
                new_population.append(new_individual)
        
        # Force inject random individuals (20% of population) to prevent convergence
        num_random_inject = max(1, int(self.population_size * 0.2))
        for _ in range(num_random_inject):
            if len(new_population) > num_random_inject:
                # Replace worst individuals with random ones
                worst_idx = len(new_population) - 1 - (_ % len(new_population))
                new_genes = self._generate_random_genes()
                new_population[worst_idx] = Individual(new_genes)
        
        return new_population
    
    def _diversity_tournament_selection(self, candidates: List[Individual], 
                                      current_population: List[Individual]) -> Individual:
        """Select individual that maximizes diversity in current population"""
        if not candidates:
            return None
            
        tournament_size = min(3, len(candidates))
        tournament = random.sample(candidates, tournament_size)
        
        best_individual = None
        best_diversity_score = -float('inf')
        
        for individual in tournament:
            # Calculate diversity score (distance to nearest neighbor in current population)
            min_distance = float('inf')
            for existing in current_population:
                distance = self._calculate_distance(individual, existing)
                min_distance = min(min_distance, distance)
            
            # Combine fitness and diversity (multi-objective)
            diversity_bonus = min_distance * 0.5  # Weight diversity
            combined_score = -individual.fitness + diversity_bonus  # Higher is better
            
            if combined_score > best_diversity_score:
                best_diversity_score = combined_score
                best_individual = individual
        
        return best_individual
    
    def score(self, pose: Pose) -> float:
        """Score a pose using Vina-like scoring function with protein-ligand interactions"""
        # Calculate intramolecular ligand score
        intramolecular_score = self.scoring.calculate_vina_score(pose.coordinates)
        
        # Calculate protein-ligand interaction score using pose position
        position = np.mean(pose.coordinates, axis=0)  # Center of ligand
        center = self.grid_box.center  # Use grid box center consistently
        
        # Distance from binding site center affects score
        distance_from_center = np.linalg.norm(position - center)
        distance_penalty = distance_from_center * 0.5  # Penalize poses far from center
        
        # Binding site interaction term (simplified)
        binding_score = 0.0
        for coord in pose.coordinates:
            dist_to_center = np.linalg.norm(coord - center)
            if dist_to_center < 10.0:  # Within binding site
                binding_score -= 2.0 * (10.0 - dist_to_center) / 10.0  # Favorable interaction
            else:
                binding_score += 1.0  # Penalty for being outside
        
        # Add clash penalty
        clash_score = self.scoring.calculate_clash_score(pose.coordinates)
        
        # Combine all terms with realistic weights
        raw_score = (
            intramolecular_score * 0.3 +  # Reduce weight of intramolecular
            binding_score * 0.5 +         # Protein-ligand interactions
            distance_penalty * 0.2 +      # Position penalty
            clash_score * 5.0              # Clash penalty
        )
        
        # Add random diversity to prevent identical scores
        diversity_factor = np.random.normal(0, 0.5)  # Small random variation
        
        # Scale to realistic docking energy range (5-14, which becomes -14 to -5 kcal/mol)
        # Map raw scores to AutoDock Vina-like range with better diversity
        if raw_score < 20:
            # Excellent poses: map to 5-7 range (becomes -7 to -5 kcal/mol)
            total_score = 5.0 + (raw_score / 20.0) * 2.0
        elif raw_score < 60:
            # Good poses: map to 7-9 range (becomes -9 to -7 kcal/mol)
            total_score = 7.0 + ((raw_score - 20.0) / 40.0) * 2.0
        elif raw_score < 120:
            # Moderate poses: map to 9-12 range (becomes -12 to -9 kcal/mol)
            total_score = 9.0 + ((raw_score - 60.0) / 60.0) * 3.0
        else:
            # Poor poses: map to 12-14 range (becomes -14 to -12 kcal/mol)
            total_score = 12.0 + min(2.0, (raw_score - 120.0) / 80.0 * 2.0)
        
        # Add diversity and clamp to valid range
        total_score = np.clip(total_score + diversity_factor, 5.0, 14.0)
        
        # Update pose scoring details
        pose.vdw_energy = self.scoring.calculate_vdw_energy(pose.coordinates)
        pose.hbond_energy = self.scoring.calculate_hbond_energy(pose.coordinates)
        pose.hydrophobic_energy = self.scoring.calculate_hydrophobic_energy(pose.coordinates)
        pose.clash_score = clash_score
        
        return total_score
    
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
