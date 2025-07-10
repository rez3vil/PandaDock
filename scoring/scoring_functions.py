"""
Comprehensive scoring functions for molecular docking
Implements various energy terms and empirical scoring functions
"""

import numpy as np
from typing import List, Dict, Any, Optional, Tuple
import logging

from ..utils.math_utils import distance_matrix, angle_between_vectors


class ScoringFunctions:
    """
    Comprehensive scoring functions for molecular docking
    
    Features:
    - Van der Waals interactions
    - Electrostatic interactions
    - Hydrogen bonding
    - Hydrophobic interactions
    - Solvation effects
    - Entropy penalties
    - Vina-style empirical scoring
    - Glide-style physics-based scoring
    """
    
    def __init__(self, config=None):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Energy parameters
        self.vdw_parameters = self._initialize_vdw_parameters()
        self.electrostatic_parameters = self._initialize_electrostatic_parameters()
        self.hbond_parameters = self._initialize_hbond_parameters()
        self.hydrophobic_parameters = self._initialize_hydrophobic_parameters()
        self.solvation_parameters = self._initialize_solvation_parameters()
        
        # Scoring weights
        self.weights = {
            'vdw': config.scoring.vdw_weight if config else 1.0,
            'electrostatic': config.scoring.electrostatic_weight if config else 1.0,
            'hbond': config.scoring.hbond_weight if config else 1.0,
            'hydrophobic': config.scoring.hydrophobic_weight if config else 1.0,
            'solvation': config.scoring.solvation_weight if config else 1.0,
            'entropy': config.scoring.entropy_weight if config else 1.0
        }
        
        # Cutoff distances
        self.vdw_cutoff = 8.0  # Angstroms
        self.electrostatic_cutoff = 12.0  # Angstroms
        self.hbond_cutoff = 3.5  # Angstroms
        self.hydrophobic_cutoff = 4.5  # Angstroms
        
        self.logger.info("Initialized ScoringFunctions")
    
    def _initialize_vdw_parameters(self) -> Dict[str, Dict[str, float]]:
        """Initialize van der Waals parameters"""
        # Lennard-Jones parameters (epsilon in kcal/mol, sigma in Angstroms)
        return {
            'H': {'epsilon': 0.0157, 'sigma': 2.5},
            'C': {'epsilon': 0.0860, 'sigma': 3.4},
            'N': {'epsilon': 0.1700, 'sigma': 3.25},
            'O': {'epsilon': 0.2100, 'sigma': 2.96},
            'P': {'epsilon': 0.2000, 'sigma': 3.74},
            'S': {'epsilon': 0.2500, 'sigma': 3.55},
            'F': {'epsilon': 0.0610, 'sigma': 2.94},
            'Cl': {'epsilon': 0.2760, 'sigma': 3.52},
            'Br': {'epsilon': 0.3890, 'sigma': 3.73},
            'I': {'epsilon': 0.5500, 'sigma': 3.96}
        }
    
    def _initialize_electrostatic_parameters(self) -> Dict[str, float]:
        """Initialize electrostatic parameters"""
        return {
            'dielectric_constant': 4.0,  # Protein dielectric
            'coulomb_constant': 332.0636  # kcal*A/mol/e^2
        }
    
    def _initialize_hbond_parameters(self) -> Dict[str, Any]:
        """Initialize hydrogen bond parameters"""
        return {
            'donors': ['N', 'O'],
            'acceptors': ['N', 'O'],
            'optimal_distance': 2.8,  # Angstroms
            'optimal_angle': 180.0,  # degrees
            'energy_scale': 5.0,  # kcal/mol
            'distance_tolerance': 0.5,  # Angstroms
            'angle_tolerance': 30.0  # degrees
        }
    
    def _initialize_hydrophobic_parameters(self) -> Dict[str, Any]:
        """Initialize hydrophobic interaction parameters"""
        return {
            'hydrophobic_atoms': ['C'],
            'optimal_distance': 3.8,  # Angstroms
            'energy_scale': 0.5,  # kcal/mol
            'distance_tolerance': 1.0  # Angstroms
        }
    
    def _initialize_solvation_parameters(self) -> Dict[str, float]:
        """Initialize solvation parameters"""
        # Atomic solvation parameters (kcal/mol/A^2)
        return {
            'H': 0.0,
            'C': 0.0118,
            'N': -0.0015,
            'O': -0.0009,
            'P': 0.0,
            'S': 0.0118,
            'F': -0.0009,
            'Cl': 0.0118,
            'Br': 0.0118,
            'I': 0.0118
        }
    
    def calculate_total_energy(self, ligand_coords: np.ndarray, protein_coords: np.ndarray = None) -> float:
        """
        Calculate total energy for a ligand pose
        
        Args:
            ligand_coords: Ligand coordinates
            protein_coords: Protein coordinates (optional)
            
        Returns:
            Total energy in kcal/mol
        """
        total_energy = 0.0
        
        # Intramolecular ligand energy
        if len(ligand_coords) > 1:
            total_energy += self.calculate_intramolecular_energy(ligand_coords)
        
        # Intermolecular ligand-protein energy
        if protein_coords is not None:
            total_energy += self.calculate_intermolecular_energy(ligand_coords, protein_coords)
        
        return total_energy
    
    def calculate_intramolecular_energy(self, coordinates: np.ndarray) -> float:
        """Calculate intramolecular energy within ligand"""
        if len(coordinates) < 2:
            return 0.0
        
        # Calculate distances
        distances = distance_matrix(coordinates, coordinates)
        
        # Remove diagonal and upper triangle
        mask = np.triu(np.ones(distances.shape), k=1).astype(bool)
        
        energy = 0.0
        
        # Van der Waals energy
        for i in range(len(coordinates)):
            for j in range(i + 1, len(coordinates)):
                if distances[i, j] < self.vdw_cutoff:
                    # Skip 1-2 and 1-3 interactions (would need bond information)
                    if j - i > 2:  # Simplified
                        vdw_energy = self._calculate_vdw_pair_energy(
                            distances[i, j], 'C', 'C'  # Simplified atom types
                        )
                        energy += vdw_energy
        
        return energy
    
    def calculate_intermolecular_energy(self, ligand_coords: np.ndarray, protein_coords: np.ndarray) -> float:
        """Calculate intermolecular energy between ligand and protein"""
        total_energy = 0.0
        
        # Calculate all pairwise distances
        distances = distance_matrix(ligand_coords, protein_coords)
        
        # Van der Waals energy
        vdw_energy = 0.0
        for i in range(len(ligand_coords)):
            for j in range(len(protein_coords)):
                dist = distances[i, j]
                if dist < self.vdw_cutoff:
                    vdw_energy += self._calculate_vdw_pair_energy(dist, 'C', 'C')  # Simplified
        
        total_energy += vdw_energy * self.weights['vdw']
        
        # Electrostatic energy
        electrostatic_energy = 0.0
        for i in range(len(ligand_coords)):
            for j in range(len(protein_coords)):
                dist = distances[i, j]
                if dist < self.electrostatic_cutoff:
                    electrostatic_energy += self._calculate_electrostatic_pair_energy(
                        dist, 0.0, 0.0  # Simplified charges
                    )
        
        total_energy += electrostatic_energy * self.weights['electrostatic']
        
        return total_energy
    
    def calculate_vdw_energy(self, coordinates: np.ndarray) -> float:
        """Calculate van der Waals energy"""
        if len(coordinates) < 2:
            return 0.0
        
        distances = distance_matrix(coordinates, coordinates)
        energy = 0.0
        
        for i in range(len(coordinates)):
            for j in range(i + 1, len(coordinates)):
                dist = distances[i, j]
                if dist < self.vdw_cutoff:
                    energy += self._calculate_vdw_pair_energy(dist, 'C', 'C')  # Simplified
        
        return energy
    
    def _calculate_vdw_pair_energy(self, distance: float, atom_type1: str, atom_type2: str) -> float:
        """Calculate van der Waals energy for atom pair"""
        if distance <= 0:
            return 1000.0  # Large penalty for zero distance
        
        # Get Lennard-Jones parameters
        params1 = self.vdw_parameters.get(atom_type1, self.vdw_parameters['C'])
        params2 = self.vdw_parameters.get(atom_type2, self.vdw_parameters['C'])
        
        # Combining rules
        epsilon = np.sqrt(params1['epsilon'] * params2['epsilon'])
        sigma = (params1['sigma'] + params2['sigma']) / 2
        
        # Lennard-Jones potential
        r_over_sigma = sigma / distance
        r6 = r_over_sigma ** 6
        r12 = r6 ** 2
        
        energy = 4 * epsilon * (r12 - r6)
        
        return energy
    
    def calculate_electrostatic_energy(self, coordinates: np.ndarray) -> float:
        """Calculate electrostatic energy"""
        # This is a simplified implementation
        # In practice, would use actual partial charges
        return 0.0
    
    def _calculate_electrostatic_pair_energy(self, distance: float, charge1: float, charge2: float) -> float:
        """Calculate electrostatic energy for atom pair"""
        if distance <= 0 or charge1 == 0 or charge2 == 0:
            return 0.0
        
        # Coulomb potential with distance-dependent dielectric
        dielectric = self.electrostatic_parameters['dielectric_constant']
        coulomb_k = self.electrostatic_parameters['coulomb_constant']
        
        energy = coulomb_k * charge1 * charge2 / (dielectric * distance)
        
        return energy
    
    def calculate_hbond_energy(self, coordinates: np.ndarray) -> float:
        """Calculate hydrogen bond energy"""
        # This is a simplified implementation
        # In practice, would identify actual H-bond donors/acceptors
        
        if len(coordinates) < 3:
            return 0.0
        
        distances = distance_matrix(coordinates, coordinates)
        energy = 0.0
        
        # Find potential H-bonds
        for i in range(len(coordinates)):
            for j in range(i + 1, len(coordinates)):
                dist = distances[i, j]
                if dist < self.hbond_cutoff:
                    # Simplified H-bond energy
                    optimal_dist = self.hbond_parameters['optimal_distance']
                    energy_scale = self.hbond_parameters['energy_scale']
                    
                    if dist < optimal_dist + 0.5:
                        energy -= energy_scale * np.exp(-(dist - optimal_dist)**2 / 0.25)
        
        return energy
    
    def calculate_hydrophobic_energy(self, coordinates: np.ndarray) -> float:
        """Calculate hydrophobic interaction energy"""
        # This is a simplified implementation
        # In practice, would identify actual hydrophobic atoms
        
        if len(coordinates) < 2:
            return 0.0
        
        distances = distance_matrix(coordinates, coordinates)
        energy = 0.0
        
        # Find potential hydrophobic interactions
        for i in range(len(coordinates)):
            for j in range(i + 1, len(coordinates)):
                dist = distances[i, j]
                if dist < self.hydrophobic_cutoff:
                    # Simplified hydrophobic energy
                    optimal_dist = self.hydrophobic_parameters['optimal_distance']
                    energy_scale = self.hydrophobic_parameters['energy_scale']
                    
                    if dist < optimal_dist + 1.0:
                        energy -= energy_scale * np.exp(-(dist - optimal_dist)**2 / 1.0)
        
        return energy
    
    def calculate_solvation_energy(self, coordinates: np.ndarray) -> float:
        """Calculate solvation energy"""
        # This is a simplified implementation
        # In practice, would calculate solvent accessible surface area
        
        if len(coordinates) < 2:
            return 0.0
        
        # Estimate buried surface area
        total_surface_area = len(coordinates) * 20.0  # Simplified
        buried_area = min(total_surface_area * 0.5, 100.0)  # Simplified
        
        # Solvation energy proportional to buried area
        solvation_energy = buried_area * 0.005  # kcal/mol per A^2
        
        return -solvation_energy  # Negative because burying surface is favorable
    
    def calculate_entropy_penalty(self, coordinates: np.ndarray) -> float:
        """Calculate entropy penalty for ligand binding"""
        # This is a simplified implementation
        # In practice, would consider rotatable bonds and conformational entropy
        
        num_atoms = len(coordinates)
        
        # Estimate number of rotatable bonds
        estimated_rotatable_bonds = max(0, num_atoms // 4 - 1)
        
        # Entropy penalty per rotatable bond
        entropy_penalty_per_bond = 0.6  # kcal/mol
        
        total_entropy_penalty = estimated_rotatable_bonds * entropy_penalty_per_bond
        
        return total_entropy_penalty
    
    def calculate_clash_score(self, coordinates: np.ndarray) -> float:
        """Calculate clash score based on van der Waals overlaps"""
        if len(coordinates) < 2:
            return 0.0
        
        distances = distance_matrix(coordinates, coordinates)
        clash_score = 0.0
        
        # Van der Waals radii
        vdw_radius = 1.7  # Simplified (carbon radius)
        
        for i in range(len(coordinates)):
            for j in range(i + 1, len(coordinates)):
                dist = distances[i, j]
                min_distance = 2 * vdw_radius
                
                if dist < min_distance:
                    overlap = min_distance - dist
                    clash_score += overlap * overlap  # Quadratic penalty
        
        return clash_score
    
    def calculate_vina_score(self, coordinates: np.ndarray) -> float:
        """
        Calculate Vina-style empirical scoring function
        
        Vina uses a knowledge-based scoring function with terms for:
        - Gauss interactions
        - Repulsive interactions
        - Hydrophobic interactions
        - Hydrogen bonds
        - Rotatable bond penalty
        """
        if len(coordinates) < 2:
            return 0.0
        
        distances = distance_matrix(coordinates, coordinates)
        
        # Vina scoring terms
        gauss1_score = 0.0
        gauss2_score = 0.0
        repulsive_score = 0.0
        hydrophobic_score = 0.0
        hbond_score = 0.0
        
        # Vina parameters
        gauss1_offset = 0.0
        gauss1_width = 0.5
        gauss2_offset = 3.0
        gauss2_width = 2.0
        
        for i in range(len(coordinates)):
            for j in range(i + 1, len(coordinates)):
                dist = distances[i, j]
                
                if dist < 8.0:  # Cutoff distance
                    # Gauss 1 term
                    gauss1_score += np.exp(-((dist - gauss1_offset) / gauss1_width)**2)
                    
                    # Gauss 2 term
                    gauss2_score += np.exp(-((dist - gauss2_offset) / gauss2_width)**2)
                    
                    # Repulsive term
                    if dist < 0.0:
                        repulsive_score += dist * dist
                    
                    # Hydrophobic term (simplified)
                    if 1.5 < dist < 4.5:
                        hydrophobic_score += 1.0
                    
                    # Hydrogen bond term (simplified)
                    if 1.0 < dist < 3.5:
                        hbond_score += 1.0
        
        # Rotatable bond penalty
        rotatable_bonds = max(0, len(coordinates) // 4 - 1)  # Simplified
        rotatable_penalty = rotatable_bonds * 0.5854
        
        # Combine terms (Vina weights)
        total_score = (
            -0.0356 * gauss1_score +
            -0.00516 * gauss2_score +
            0.840 * repulsive_score +
            -0.0351 * hydrophobic_score +
            -0.587 * hbond_score +
            rotatable_penalty
        )
        
        return total_score
    
    def find_hbond_interactions(self, coordinates: np.ndarray) -> List[Dict[str, Any]]:
        """Find hydrogen bond interactions"""
        interactions = []
        
        if len(coordinates) < 3:
            return interactions
        
        distances = distance_matrix(coordinates, coordinates)
        
        # Find potential H-bonds
        for i in range(len(coordinates)):
            for j in range(i + 1, len(coordinates)):
                dist = distances[i, j]
                if dist < self.hbond_cutoff:
                    # Simplified H-bond detection
                    interaction = {
                        'donor_atom': i,
                        'acceptor_atom': j,
                        'distance': dist,
                        'angle': 180.0,  # Simplified
                        'energy': -2.0 * np.exp(-(dist - 2.8)**2 / 0.25)
                    }
                    interactions.append(interaction)
        
        return interactions
    
    def find_hydrophobic_interactions(self, coordinates: np.ndarray) -> List[Dict[str, Any]]:
        """Find hydrophobic interactions"""
        interactions = []
        
        if len(coordinates) < 2:
            return interactions
        
        distances = distance_matrix(coordinates, coordinates)
        
        # Find potential hydrophobic interactions
        for i in range(len(coordinates)):
            for j in range(i + 1, len(coordinates)):
                dist = distances[i, j]
                if dist < self.hydrophobic_cutoff:
                    # Simplified hydrophobic detection
                    interaction = {
                        'atom1': i,
                        'atom2': j,
                        'distance': dist,
                        'energy': -0.5 * np.exp(-(dist - 3.8)**2 / 1.0)
                    }
                    interactions.append(interaction)
        
        return interactions
    
    def find_salt_bridge_interactions(self, coordinates: np.ndarray) -> List[Dict[str, Any]]:
        """Find salt bridge interactions"""
        interactions = []
        
        # This is a simplified implementation
        # In practice, would identify actual charged atoms
        
        if len(coordinates) < 2:
            return interactions
        
        distances = distance_matrix(coordinates, coordinates)
        
        # Find potential salt bridges
        for i in range(len(coordinates)):
            for j in range(i + 1, len(coordinates)):
                dist = distances[i, j]
                if dist < 4.0:  # Salt bridge cutoff
                    # Simplified salt bridge detection
                    interaction = {
                        'positive_atom': i,
                        'negative_atom': j,
                        'distance': dist,
                        'energy': -5.0 * np.exp(-(dist - 3.0)**2 / 0.5)
                    }
                    interactions.append(interaction)
        
        return interactions
    
    def get_scoring_info(self) -> Dict[str, Any]:
        """Get information about scoring functions"""
        return {
            'scoring_terms': [
                'van_der_waals',
                'electrostatic',
                'hydrogen_bonding',
                'hydrophobic',
                'solvation',
                'entropy'
            ],
            'weights': self.weights,
            'cutoff_distances': {
                'vdw': self.vdw_cutoff,
                'electrostatic': self.electrostatic_cutoff,
                'hbond': self.hbond_cutoff,
                'hydrophobic': self.hydrophobic_cutoff
            },
            'vina_terms': [
                'gauss1',
                'gauss2',
                'repulsive',
                'hydrophobic',
                'hydrogen_bond',
                'rotatable_bond_penalty'
            ]
        }